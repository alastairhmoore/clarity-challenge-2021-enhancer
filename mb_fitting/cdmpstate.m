function [next,ixu]=cdmpstate(n)
% create quadratic programming state machine
%  Inputs: n    number of fixed basis elements
%
% Outputs:  next(k,2^n)  next state where 2nd argument is 1+pow2(n-1:-1:0)*(s<0)
%                        next=0 when no tests remaining, next=-1 for an error
%                           # free variables = n: 0 1 2  3  4
%                                   # states = k: 0 1 4 25 678
%           ixu(k,1)     index into combination of unconstrained variables in range 1:2^n
%                        ixu(*)-1 is a bitmap corresponding to msku(*,:)
%
% Note: all length-n bitmaps have the MSB corresponding to s(1)
persistent n0 next0 ixu0
if n==n0
    next=next0;
    ixu=ixu0;
else
    if n==0
        next=[];
        msku=[];
        ixu=zeros(0,1);
    else
        m=pow2(n); % Number of permutations of constrained and unconstrained variables
        opts=(pow2(n)-1:-1:1)'; % list option bitmaps: 1=unconstrained, 0=constrained to zero
        bmweights=pow2(1-n:0);
        bmall=m-1; % bitmap with all variables selected
        [con,ix]=sort(sum(rem(floor(opts*bmweights),2)==0,2)); % sort by number of constrained variables
        opts=opts(ix); % sorted into decreasing number of constrained variables
        k=1; % one state initially
        mskf=true(1,m-1); % all combinations are possible initially
        ixu=zeros(1,1); % bitmap of unconstrained variables to evaluate in each state
        msku=true(1,n); % same as bmu but expanded to a mask
        next=zeros(1,m); % next state for each possibile value of 1+pow2(n-1:-1:0)*(s<0)
        i=1; % current state
        while i<=k
            itrial=find(mskf(i,:),1); % find the highest priority combination in this state
            bmui=opts(itrial); % bitmap of trial unconstrained variables
            ixu(i,1)=bmui+1;
            msku(i,:)=rem(floor(bmui*bmweights),2); % convert to a mask also
            for j=1:m                           % for every possible combination of (s<0)
                mskfj=mskf(i,:);                % feasible combinatins in the current state
                bmslt=j-1;                      % bitmap of (s<0)
                bmsge=bmall-bmslt;              % bitmap of (s>=0)
                if bitor(bmui,bmsge)~=bmall     % some constrained variables are <0: this situation cannot happen
                    next(i,j)=-1;               % -1 indicates error if we reach this branch of the state machine
                else
                    mskfj(itrial)=false;        % don't test this combination again
                    if bmslt~=0             % if there are any violating variables
                        % bitand(bmall-bmui,opts)==0 is true for options that constrain all variables that are constrained in the current state
                        % bmslt identifies violating variables
                        % bitand(bmslt,bmall-opts) identifies violating variables that are constrained in an option
                        % bitand(bmslt,bmall-opts)==0 is true for options in which at none of the violating variables is constrained
                        mskfj(bitand(bmall-bmui,opts)==0 &  bitand(bmslt,bmall-opts)==0)=false; % at least one violating variable must be constrained to zero
                    end
                    if ~any(mskfj)
                        next(i,j)=0;            % no combinations left to test
                    else
                        jx=find(all(mskf==repmat(mskfj,k,1),2),1); % find if an existing state matches the remaiing possibilities
                        if isempty(jx)          % else create a new state
                            k=k+1;              % increment number of states
                            mskf(k,:)=mskfj;    % save the new list of feasible constraint combinations
                            next(i,j)=k;        % mark the next state as this new one
                        else                    % go to an existing state
                            next(i,j)=jx;
                        end
                    end
                end
            end
            i=i+1; % now process the next state
        end % end while loop
    end
    n0=n; % save in cache for next call with the same n
    next0=next;
    ixu0=ixu;
end
if ~nargout
    fprintf('State Free   Bitmap of s<0\n     Bitmap');
    fprintf('%3d ',0:size(next,2)-1);
    fprintf('\n-----------\n');
    for i=1:length(ixu)
        fprintf(' %3d %4d ',i,ixu(i)-1);
        fprintf(' %3d',next(i,:));
        fprintf('\n');
    end
end
