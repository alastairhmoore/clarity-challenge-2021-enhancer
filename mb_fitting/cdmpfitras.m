function [rmft,sdft,doft,erft]=cdmpfitras(ry,c,vq)
% Fit a multicomponent model to a microphone-array covariance matrix
%  Inputs: ry(nm^2,nt,nf)   Smoothed real-vectorized observed covariance matrices
%          c(nm^2,nb-1,nf)  Real-vectorized covariance matrices for fixed components
%          vq(nm^2,nq,nf)   Real-vectorized covariance matrices for directional component [omit or use [] if unwanted
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
%           case when the Lagrange multiplier is slightly negative due to rounding errors (but this
%           modification has not yet been added into the diretion-independent search).
% 2020_0619 Includes the output rmft() which gives the estimated version of ry.
% 2020_0711 Solves the quadratic programming problem using the algorithm from [2]. It is applied to
%           the dual problem to find the direction-independent solution and then, using this as a
%           starting point, to the primal problem to include a directional component.
% 2020_0711 Speeded up by pre-calculating Hiq
% 2020_0717 Added code to allow vq=[] if no directional component
% 2020_0731 Added code to cope with duplicate columns in c and vq
% 2020_0801 Fixed error when vq input argument is missing
% 2020_0803 Fixed error when some but not all columns of c are replicated
%
% Refs:
% [1] P. A. Naylor, A. H. Moore, and M. Brookes.
%     Improving robustness of adaptive beamforming for hearing devices.
%     In Proc Intl Symp on Auditory and Audiological Research, Nyborg, Denmark, Aug. 2019.
% [2] P. E. Gill and E. Wong. Methods for convex and general quadratic programming.
%     Mathematical Programming Computation, 7 (1): 71–112, Aug 2014. doi: 10.1007/s12532-014-0075-x.
%
[nm2,nt,nf]=size(ry);       % [#microphones^2, #time-frames, #frequency-bins]
nbf=size(c,2);              % # fixed basis-components
nb=nbf+1;                   % # total basis-components
%
rankth=1e-6;                % threshold for determining rank
doft=zeros(nt,nf);          % optimal direction per TF bin
erft=zeros(nt,nf);          % squared error per TF bin
if nargin<3 || isempty(vq)      % no directional components to find
    sdft=zeros(nbf,nt,nf);     	% optimal weights per TF bin
    for jf=1:nf                	% loop for each frequency bin
        cjf=c(:,:,jf);        	% fixed basis elements
        [qcjf,rcjf,ecjf]=qr(cjf,0); % QR decomposition of cjf to detemine its rank
        rcjfd=abs(diag(rcjf));
        rcjfz=rcjfd>rankth*rcjfd(1); % detect near-zero diagonal elements in R matrix
        cjf=cjf(:,ecjf(rcjfz)); % select an independent set of columns from cjf
        nbfz=size(cjf,2); % number of independent columns in cjf
        nbz=nbfz+1;
        Hu=rcjf(1:nbfz,1:nbfz)'*rcjf(1:nbfz,1:nbfz); % = cjf'*cjf but more efficient
        Hui=inv(Hu);
        su0=zeros(nbfz,1);           % default to all zero weights
         eeu=eye(nbfz);             	% useful identity matrix
        for it=1:nt             % loop for each frame (could maybe vectorize this to speed it up)
            %         v_finishat([jf nf; it-1 nt]); % uncomment to print out how long it will take
            ryft=ry(:,it,jf);   % smoothed observation for this T-F cell
            tu=cjf'*ryft;
            ti=Hu\tu;
            % iteration 0
            n=false(nbfz,1); 	% mask of constrained dimensinons
            mm=su0;             % initial optimal point: all Lagrange multipliers zero
            su=ti;              % initial optimal point: all weights from Hu\t
            [xs,s]=min(su);   	% find most negative Lagrange multiplier
            if xs<0             % quit if all Lagrange multipliers are non-negative
                % iteration 1
                mm(s)=-ti(s)/Hui(s,s);	% find unconstrained minimum of dimension s
                n(s)=true;              % dimension s is now constrained
                su=Hui(:,s)*mm(s)+ti;
                su(s)=0;               	% force basic set multiplier to zero
                [xs,s]=min(su);        	% find most negative Lagrange multiplier
                while xs<0              % loop while any Lagrange multipliers are negative
                    while true          % loop until m(s) is driven to zero
                        p=eeu(:,s);     % initialize to the appropriate column of the identity matrix
                        p(n)=Hui(n,n)\-Hui(n,s);    % find direction
                        ap=-xs/(Hui(s,:)*p);        % distance to minimum
                        mp=mm+ap*p;                 % try going all the way to the minimum
                        h=p<0 & mp<0 & mm>0;        % possible problem dimensions
                        if any(h)
                            ih=find(h);             % list of problem dimensions
                            [vx,sx]=min(mp(h)./mm(h));  % limiting dimension is ih(sx)
                            jh=ih(sx);              % index of limiting variable
                            mm=mm+ap/(1-vx)*p;      % update mm
                            xs=Hui(s,:)*mm+ti(s);  	% update Lagrange multipliers
                            n(jh)=false;            % mark limiting dimension as unconstrained in m
                            mm(jh)=0;               % and force mm component to exactly zero
                        else
                            mm=mp;
                            n(s)=true;              % dimension s is now constrained
                            break;
                        end
                    end
                    su(~n)=Hui(~n,n)*mm(n)+ti(~n);
                    su(n)=0;                            % force basic set multipliers to zero
                    [xs,s]=min(su);                     % find most negative Lagrange multiplier
                end
            end
            sdft(ecjf(rcjfz),it,jf)=su;                 % optimal weights for selected components
            doft(it,jf)=0;                              % optimal direction index
            rmftij=cjf*su;                              % model output
            rmft(:,it,jf)= rmftij;                      % save model output
            err=rmftij-ryft;                            % fitting error
            erft(it,jf)=err'*err;                       % squared Frobenius norm
        end
    end
else
    nq=size(vq,2);              % #search-directions
    sdft=zeros(nb,nt,nf);       % optimal weights per TF bin
    for jf=1:nf                	% loop for each frequency bin
        vqjf=vq(:,:,jf);       	% vectorized directional covariance matrix [nm2 nq]
        cjf=c(:,:,jf);        	% fixed basis elements
        [qcjf,rcjf,ecjf]=qr(cjf,0); % QR decomposition of cjf to detemine its rank
        rcjfd=abs(diag(rcjf));
        rcjfz=rcjfd>rankth*rcjfd(1); % detect near-zero diagonal elements in R matrix
        cjf=cjf(:,ecjf(rcjfz)); % select an independent set of columns from cjf
        nbfz=size(cjf,2); % number of independent columns in cjf
        nbz=nbfz+1;
        Hu=rcjf(1:nbfz,1:nbfz)'*rcjf(1:nbfz,1:nbfz); % = cjf'*cjf but more efficient
        Hui=inv(Hu);
        cvqjf=cjf'*vqjf;
        Hiq=zeros(nbz,nbz,nq);        % space for quadratic-term matrices for primal problem
        Hiq(2:nbz,2:nbz,:)=repmat(Hu,[1 1 nq]);   % fill in entries of
        Hiq(1,1,:)=sum(vqjf.^2,1);             	% ... mx=vqjf(:,iq); Hiq(:,:,iq)=[mx cjf]'*[mx cjf];
        Hiq(2:nbz,1,:)=cvqjf;
        Hiq(1,2:nbz,:)=cvqjf;
        dchendvq=Hu\cvqjf;
        kd=sum(vqjf.^2,1);                      % directional squared Frobenius norms [1 nq]
        detdi=1./(kd-sum(cvqjf.*dchendvq,1)); % normalizing determinant reciprocal [1 nq]
        pd=(vqjf-cjf*dchendvq).*repmat(sqrt(detdi),nm2,1); % calculate p vector for each direction [nm2 nq]
        su0=zeros(nbfz,1);         	% default to all zero weights
        eeu=eye(nbfz);             	% useful identity matrix
        ee=eye(nbz);               	% useful identity matrix for primal problem
        m=zeros(nbz,1);           	% default to all-zero Lagrange multipliers for primal problem
        t=zeros(nbz,1);             % inner product of basis and target vectors
        for it=1:nt                 % loop for each frame (could maybe vectorize this to speed it up)
            %         v_finishat([jf nf; it-1 nt]); % uncomment to print out how long it will take
            ryft=ry(:,it,jf);   % smoothed observation for this T-F cell
            tu=cjf'*ryft; % inner product of basis and target vectors
            ti=Hu\tu;
            %
            % First find the best solution without any directional component
            % by solving the dual problem with all dimensions initially unconstrained
            %
            % iteration 0
            n=false(nbfz,1);     % mask of constrained dimensinons
            mm=su0;             % initial optimal point: all Lagrange multipliers zero
            su=ti;              % initial optimal point: all weights from Hu\t
            [xs,s]=min(su);   	% find most negative x-component
            if xs<0             % quit if all x-components are non-negative
                % iteration 1
                mm(s)=-ti(s)/Hui(s,s);	% find unconstrained minimum of dimension s
                n(s)=true;              % dimension s is now constrained
                su=Hui(:,s)*mm(s)+ti;
                su(s)=0;               	% force basic set multiplier to zero
                [xs,s]=min(su);        	% find most negative x-component
                while xs<0              % loop while any x-components are negative
                    while true          % loop until Lagrange multiplier m(s) is driven to zero
                        p=eeu(:,s);     % initialize to the appropriate column of the identity matrix
                        p(n)=Hui(n,n)\-Hui(n,s);    % find direction
                        ap=-xs/(Hui(s,:)*p);        % distance to minimum
                        mp=mm+ap*p;                 % try going all the way to the minimum
                        h=p<0 & mp<0 & mm>0;        % possible problem dimensions
                        if any(h)
                            ih=find(h);             % list of problem dimensions
                            [vx,sx]=min(mp(h)./mm(h));  % limiting dimension is ih(sx)
                            jh=ih(sx);              % index of limiting variable
                            mm=mm+ap/(1-vx)*p;      % update Lagrange multipliers
                            xs=Hui(s,:)*mm+ti(s);  	% update x-components
                            n(jh)=false;            % mark limiting dimension as unconstrained in m
                            mm(jh)=0;               % and force Lagrange multiplier to exactly zero
                        else
                            mm=mp;
                            n(s)=true;              % dimension s is now constrained
                            break;
                        end
                    end
                    su(~n)=Hui(~n,n)*mm(n)+ti(~n);
                    su(n)=0;                         % force basic set x-components to zero
                    [xs,s]=min(su);                  % find most negative x-component
                end
            end
            eru=(cjf*su-ryft);         	% calculate error vector
            %         fvu=eru'*eru;   	% this is the lowest error without a directional component (not needed except when debugging)
            lu=vqjf'*eru;               % Lagrange multipliers for directional components: only need consider directional component iq when lu(iq)<0
            %
            % Now add in a directional component by solvng the primal problem
            % initialized by the non-directional solution
            %
            sd0=[0; su];                    % full list of weights with zero directional component
            b0=[false; ~n];                	% initial set of unconstrained dimensions
            mk=find(lu<0 & abs(kd.*detdi*rankth)'<1);  % list of feasible directions (excl ones with singular Hiq)
            if ~isempty(mk)                 % if any feasible directions
                [v,dopti]=max(abs(ryft'*pd(:,mk))); % find optimal direction assuming all components are used
                dopt=mk(dopti);             % and the optimal direction index
                voc=[vqjf(:,dopt) cjf];
                sd=(voc'*voc)\(voc'*ryft); 	% and the weights for directional solution using all components
                t(2:nbz)=tu;
                if any(sd<0)                % KT conditions not OK for optimal direction, so some of the fixed components must be zero
                    nmk=length(mk);         % number of feasible directions
                    sdq=zeros(nbz,nq);       % weights for each direction
                    fvq=zeros(1,nq);        % fitting error for each direction
                    for im=1:nmk            % loop through all feasible directions
                        iq=mk(im);                      % index of this search direction
                        mx=vqjf(:,iq);               	% directional basis vector for this direction
                        H=Hiq(:,:,iq);                  % quadratic term matrix for augemnted basis set H=[mx cjf]'*[mx cjf];
                        t(1)=mx'*ryft;                  % and linear term vector
                        x=sd0;                         	% initialize with direction independent solution
                        b=b0;                         	% initial set of unconstrained dimensions
                        s=1;                            % investigate dimension 1
                        ms=lu(iq);                      % same as calculating directly from m=H(1,:)*x-t(1)
                        while ms<0                      % loop while any Lagrange multipliers are negative
                            while true                  % loop until m(s) is driven to zero
                                p=ee(:,s);              % initialize to the appropriate column of the identity matrix
                                p(b)=H(b,b)\-H(b,s);    % find direction
                                q=H*p;                  % change in Lagrange multipliers due to p
                                ap=-ms/(H(s,:)*p);      % distance to minimum
                                xp=x+ap*p;              % try going all the way to the minimum
                                h=p<0 & xp<0 & x>0;   	% possible problem dimensions
                                if any(h)               % check if any dimensions have become negative
                                    ih=find(h);
                                    [vx,sx]=min(xp(h)./x(h));
                                    x=x+ap/(1-vx)*p;    % update x
                                    jh=ih(sx);          % index of limiting variable
                                    b(jh)=false;        % mark limiting dimension as constrained
                                    x(jh)=0;            % and force x component to exactly zero
                                    ms=H(s,:)*x-t(s);  	% update Lagrange multiplier for dimension s
                                else
                                    x=xp;
                                    b(s)=true;          % mark dimension s as unconstrained
                                    break;              % m(s) successfully driven to zero
                                end
                            end
                            m(~b)=H(~b,b)*x(b)-t(~b);   % Lagrange multipliers
                            m(b)=0;                     % force basic set multipliers to zero
                            [ms,s]=min(m);            	% find most negative Lagrange multiplier
                        end
                        sdq(:,iq)=x;                    % save the optimum weights
                        fvq(iq)=sum(([mx cjf]*x-ryft).^2);	% squared Frobenius norm of error; should always be less than fvu
                        % optionally calculate lambda which should be positive: lamq=[vqjf(:,iq) cjf]'*err
                        % compare: [sdqp,fvqp0]=quadprog([vqjf(:,iq) cjf]'*[vqjf(:,iq) cjf],-[vqjf(:,iq) cjf]'*ryft,-eye(nbz),zeros(nbz,1),[],[],[],[],[],optimset('display','off')); fvqp=2*fvqp0+ryft'*ryft; errdb=db(fvq(iq)/fvqp)/2
                    end
                    [v,dopti]=min(fvq(mk));             % find the best direction as that with lowest error
                    dopt=mk(dopti);                     % extract optimal direction
                    sd=sdq(:,dopt);                     % and weights
                end
            else                                        % no feasible directions, so ...
                sd=sd0;                                 % use the direction-independent version
                dopt=1;                                 % default to first direction
            end
            sdfull=zeros(nb,1);
            sdfull(1)=sd(1);                            % insert directional component
            sdfull(1+ecjf(rcjfz))=sd(2:end);            % insert selected fixed components
            sdft(:,it,jf)=sdfull;                      	% optimal weights
            doft(it,jf)=dopt;                           % optimal direction index
            rmftij=[vqjf(:,dopt) cjf]*sd;             	% model output
            rmft(:,it,jf)= rmftij;                      % save model output
            err=rmftij-ryft;                            % fitting error
            erft(it,jf)=err'*err;                       % squared Frobenius norm
        end
    end
end