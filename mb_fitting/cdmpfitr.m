function [rmft,sdft,doft,erft]=cdmpfitr(ry,c,vq)
% Fit a multicomponent model to a microphone-array covariance matrix
%  Inputs: ry(nm^2,nt,nf)   Smoothed real-vectorized observed covariance matrices
%          c(nm^2,nb-1,nf)  Real-vectorized covariance matrices for fixed components
%          vq(nm^2,nq,nf)   Real-vectorized covariance matrices for directional component
%
% Outputs: rmft(nm^2,nt,nf) Model fit to ry
%          sdft(nb,nt,nf)   Component weights: directional component first followed by fixed components
%          doft(nt,nf)      Index of directional component in vq array
%          erft(nt,nf)      Squared Frobenius norm of fitting error
%
%    where nm = number of microphone channels
%          nt = number of time frames
%          nf = number of frequency bins
%          nb = number of components in the model
%          nq = number of search directions
%
% This implements the fitting stage of the covariance matrix model described in [1].
% The input arguments are all in the form of real-vectorized covariance matrices; given a
% complex nm x nm Hermitian matrix, X, the real-vectorized form is Y(:) where
%             { sqrt(2)*real(X(i,j))   for i<j
%    Y(i,j) = { X(i,j)                 for i=j
%             { sqrt(2)*imag(X(i,j))   for i>j
% The functions cdmpcovr() and cdmpcovs() create real-vectorized vectors
% from covariance matrices and STFT-domain signals respectively.
%
% Versions:
%
% 2020_0417 Now supports an arbitrary number of columns of c. The directional search copes with the
% case when the Lagrange multiplier is slightly negative due to rounding errors (but this
% modification has not yet been added into the diretion-independent search).
% 2020_0619 Includes the output rmft() which gives the estimated version of ry.
% Refs:
% [1] P. A. Naylor, A. H. Moore, and M. Brookes.
%     Improving robustness of adaptive beamforming for hearing devices.
%     In Proc Intl Symp on Auditory and Audiological Research, Nyborg, Denmark, Aug. 2019.
%
[nm2,nt,nf]=size(ry);       % [#microphones^2, #time-frames, #frequency-bins]
nbf=size(c,2);              % # fixed basis-components
nb=nbf+1;                   % # total basis-components
[next,ixu]=cdmpstate(nbf);  % create state machine to search constaraint combinations
mku2bm=pow2(nbf-1:-1:0);    % conversion from mask to bitmap
mkx2bm=[0 mku2bm];          % conversion from mask to bitmap including 0 for directional component
%
nq=size(vq,2);              % #search-directions
sdft=zeros(nb,nt,nf);       % optimal weights per TF bin
doft=zeros(nt,nf);          % optimal direction per TF bin
erft=zeros(nt,nf);          % squared error per TF bin
sdx=zeros(nb,1);            % space for weights including directional component
sdb=sdx;                    % best set of weights
su=zeros(nbf,1); % space for weights exclusding directional component
su0=zeros(nbf,1); % default to all zero weights
for jf=1:nf                	% loop for each frequency bin
    vqi=vq(:,:,jf);         % vectorized directional covariance matrix [nm2 nq]
    ci=c(:,:,jf);           % fixed basis elements
    [msku,dch,qph]=cdmpfix(ci,vqi);
    mskc=~msku; % constrained variables mask
    mskx=[true(size(msku,1),1) msku]; % extended unconstrained mask (including directinal comonent)
    kd=sum(vqi.^2,1);       % directional squared Frobenius norms [1 nq]
    bd=ci'*vqi;             % product needed below [nb-1 nq]
    detdi=1./(kd-sum(bd.*(dch{end}*vqi),1)); % normalizing determinant reciprocal [1 nq]
    pd=(vqi-ci*dch{end}*vqi).*repmat(sqrt(detdi),nm2,1); % calculate p vector for each direction [nm2 nq]
    for it=1:nt             % loop for each frame (could maybe vectorize this to speed it up)
        %         v_finishat([jf nf; it-1 nt]); % uncomment to print out how long it will take
        ryft=ry(:,it,jf);   % smoothed observation for this T-F cell
        %
        % First find the best solution without any directional component
        %
        % ====== check what happens for all-constrained state ====
        is=1; % initial state
        sub=su0;  % best weights so far
        fvb=ryft'*ryft; % Squared Frobenius norm error of best so far
        while is>0
            ixui=ixu(is);                       % index of free variables (1 + bitmap)
            su(msku(ixui,:))=dch{ixui}*ryft;    % calculate constrained optimum (all elements are set when is=1)
            six=1+mku2bm*(su<0);                % 1 + bitmap of violating components
            if six==1                           % all constrains satisfied => check lagrange multipliers
                eru=(ci*su-ryft);               % error vector
                lsu=ci(:,mskc(ixui,:))'*eru;    % lagrange multipliers
                if all(lsu>=0)          % found optimum if all lsu elements >=0 or if lsu is empty (for is=1)
                    sub=su;             % save as the definitely the best solution
                    break;              % quit since we have found the minimum
                else % possibly the solution (if rounding errors make the sign of lsu wrong)
                    fvx=eru'*eru;       % find the squared Frobenius norm error
                    if fvx<fvb
                        sub=su;         % save this solution as best so far
                        fvb=fvx;        % and its squared Frobenius norm error
                    end
                end
            end
            su=su0; % set all weights to zero by default
            is=next(is,six); % next state
        end % while
        su=sub;                 % save the best solution that has been found
        eru=(ci*su-ryft);       % calculate error vector
        %         fvu=eru'*eru;           % this is the lowest error without a directional component (not needed except when debugging)
        lu=vqi'*eru;            % Lagrange multipliers for directional components: only need consider directional component iq when lu(iq)<0
        %
        % Now add in a directional component
        %
        mk=find(lu<0);                  % list of feasible directions
        if ~isempty(mk)                 % if any feasible directions
            [v,dopti]=max(abs(ryft'*pd(:,mk))); % find optimal direction assuming all components are used
            dopt=mk(dopti);             % and the optimal direction index
            sd=qph{end,dopt}*ryft;   % and the weights for directional solution using all components
            if any(sd<0)                % KT conditions not OK for optimal direction, so some of the fixed components must be zero
                nmk=length(mk);         % number of feasible directions
                sdq=zeros(nb,nq);       % weights for each direction
                fvq=zeros(1,nq);        % fitting error for each direction
                for im=1:nmk            % loop through all feasible directions
                    iq=mk(im);                      % index of this search direction
                    mx=vqi(:,iq);                   % basis vector for this direction
                    is=1; % initial state
                    sdb(:)=0; % initialize weights to zero
                    sdb(1)=max(qph{1,iq}*ryft,0); % default to all zero weights for fixed components
                    erb=mx*sdb(1)-ryft; % error vector if all fixed components constrained
                    fvb=erb'*erb; % lowest error so far
                    while is>0 % loop through state machine trying constraint combinations
                        ixui=ixu(is); % index of free variables
                        sdx(:)=0; % set all constrained weights to zero
                        sdx(mskx(ixui,:))=qph{ixui,iq}*ryft; % calculate constrained optimum
                        six=1+mkx2bm*(sdx<0); % 1 + bitmap of violating components (excluding directional component)
                        if six==1 % all constraints satisfied => check lagrange multipliers
                            erx=[mx ci]*sdx-ryft;    % calculate the fitting error vector
                            lsx=ci(:,mskc(ixui,:))'*erx; % should always be >=0 at true solution but rounding errors may prevent this
                            if all(lsx>=0) % definitely the solution (also if lsx is empty)
                                sdb=sdx; % save as the best solution
                                fvb=erx'*erx; % save squared error as well
                                break; % quit if we have definitely found the minimum
                            else % possibly the solution (if rounding errors make the sign of lsx wrong)
                                fvx=erx'*erx; % squared error
                                if fvx<fvb % check if this is the best solution so far
                                    sdb=sdx; % save this solution
                                    fvb=fvx;
                                end % now carry on checking
                            end
                        end
                        is=next(is,six); % next state
                    end % while
                    sdq(:,iq)=sdb;                  % save the optimum weights
                    fvq(iq)=fvb;               % squared Frobenius norm of error; should always be less than fvu
                    % optionally calculate lambda which should be positive: lamq=[vqi(:,iq) ci]'*err
                    % compare: [sdqp,fvqp0]=quadprog([vqi(:,iq) ci]'*[vqi(:,iq) ci],-[vqi(:,iq) ci]'*ryft,-eye(nb),zeros(nb,1),[],[],[],[],[],optimset('display','off')); fvqp=2*fvqp0+ryft'*ryft; errdb=db(fvq(iq)/fvqp)/2
                end
                [v,dopti]=min(fvq(mk));             % find the best direction as tat with lowest error
                dopt=mk(dopti);                     % extract optimal direction
                sd=sdq(:,dopt);                     % and weights
            end
        else                                        % no feasible directions, so ...
            sd=[0; su];                             % use the direction-independent version
            dopt=1;                                 % default to first direction
        end
        sdft(:,it,jf)=sd;                           % optimal weights
        doft(it,jf)=dopt;                           % optimal direction index
        rmftij=[vqi(:,dopt) ci]*sd;                 % model output
        rmft(:,it,jf)= rmftij;                      % save model output
        err=rmftij-ryft;                            % fitting error
        erft(it,jf)=err'*err;                       % squared Frobenius norm
    end
end