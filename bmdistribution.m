% Author  : Nicolas Granger <nicolas.granger@telecom-sudparis.eu>
% License : MIT (see LICENSE file)
% Source  : https://github.com/pixelou/bmdistribution

classdef bmdistribution < handle 
%BMDISTRIBUTION Binary (Bernoulli) mixture distribution class.
%   An object of the BMDISTRIBUTION class defines a Binary mixture
%   distribution, which is a multivariate distribution that consists of a
%   mixture of one or more multivariate binary (Bernoulli) distribution
%   components. The number of components for a given BMDISTRIBUTION object
%   is fixed. Each multivariate Binary component is defined by its mean,
%   and the mixture is defined by a vector of mixing proportions.
%
%   To create a binary mixture distribution by specifying the distribution
%   parameters, use the BMDISTRIBUTION constructor. To fit a binary mixture
%   distribution model to data, use FITBMDIST.
%   
%   See also FITBMDIST.

    properties(GetAccess=public, SetAccess=protected)
        %NumVariables Number of variables(dimensions).
        NumVariables = 0;
        %DistributionName Name of the distribution.
        DistributionName = 'bernoulli mixture distribution';
        %NumComponents Number of mixture components.
        NumComponents = 0;
        %ComponentProportion Mixing proportion of each component.
        ComponentProportion = zeros(0,1);
        %Means Matrix of component means.
        Means = zeros(0,1);
        %NumIterations Number of iterations.
        NumIterations = 0;
    end
    
    methods
        
        function obj = bmdistribution(M, p)
        %BMDISTRIBUTION Create a Binary (Bernoulli) mixture model.
        %   BM = BMDISTRIBUTION(M, P) creates a distribution consisting of
        %   a mixture of multivariate binary (Bernoulli) components, given
        %   values of the components' distribution parameters M and and the
        %   mixture proportions P. To create a Binary mixture distribution
        %   by fitting to data, use FITBMDIST.
        %
        %   The number of components and the dimension of the distribution are
        %   implicitly defined by the sizes of the inputs M and P.
        %
        %   M is a K-by-D matrix specifying the mean of each component,
        %   where K is the number of components, and D is the number of 
        %   dimensions.
        %
        %   P is K-wide vector specifying the mixing proportions of each 
        %   component.  If P does not sum to 1, GMDISTRIBUTION normalizes it. 
        %   The default is equal proportions if P is not given.
            assert(size(M, 1) == numel(p), ...
                'Means and priors dimensions mismatch.');
            obj.NumVariables = size(M, 2);
            obj.NumComponents = size(M, 1);
            if exist('p','var')
                p = reshape(p, numel(p), 1);
                obj.ComponentProportion = p / sum(p);
            else
                obj.ComponentProportion = ones(obj.NumComponents, 1) / obj.NumComponents;
            end
            obj.Means = M;
        end
        
        function P = posterior(this, X)
        %POSTERIOR(obj, X) Posterior probability of components given data.
        %	POST = posterior(OBJ,X) returns POST, a matrix containing
        %   estimates of the posterior probability of the components in
        %   bmdistribution OBJ given points in X. X is an N-by-D data
        %   matrix. Rows of X correspond to points, columns correspond to
        %   variables. POST(I,J) is the posterior probability of point I
        %   belonging to component J, i.e., Pr{Component J | point I}.
            K = this.NumComponents;
            N = size(X, 1);

            P = zeros(N,K);
            for k = 1:K
                P(:,k) = prod(...
                    bsxfun(@times, this.Means(k,:), X) ...
                    + bsxfun(@times, 1-this.Means(k,:), ~X), 2);
            end
        end
        
        function P = pdf(this, X)
        %PDF(OBJ, X) for a Binary (Bernoulli) mixture distribution.
        %Y = pdf(OBJ,X) returns Y, a vector of length N containing the
        %   probability density function (pdf) for the gmdistribution OBJ,
        %   evaluated at the N-by-D data matrix X. Rows of X correspond to
        %   points, columns correspond to variables. Y(I) is the pdf value
        %   of point I.
            P = this.posterior(X) * this.ComponentProportion;
        end
        
        function idx = cluster(this, X)
        %CLUSTER(OBJ, X) Cluster data for a Binary (Bernoulli) mixture
        %distribution.
        %   IDX = CLUSTER(OBJ, X) partitions the points in the N-by-D data
        %   matrix X into K clusters determined by the K components of the
        %   binary mixture distribution defined by OBJ. In the matrix X,
        %   Rows of X correspond to points, columns correspond to
        %   variables. cluster returns an N-by-1 vector IDX containing the
        %   cluster index of each point. The cluster index refers to the
        %   component giving the largest posterior probability for the
        %   point.
            P = bsxfun(@times, this.posterior(X), this.ComponentProportion');
            [~, idx] = max(P, [], 2);
            idx = reshape(idx, length(idx), 1);
        end
        
    end
    
    methods(Static=true, Hidden=true)
        
        function obj = fit(X, K, varargin)
            %FIT(X, K) adjust distribution to samples.
            %Deprectated: use fitbmdist instead
            startMethod = 'random';
            margin      = 0.1;
            maxIt       = 500;
            thres       = 10^-6;
            
            varargin = varargin{:};
            if mod(length(varargin), 2) ~= 0
                error('Incomplete Option/value pair');
            end
            i = 1;
            while i <= length(varargin)
                switch varargin{i}
                    case {'Margin'}
                        if isnumeric(varargin{i+1}) ...
                            && margin >= 0 && margin < 0.5
                            margin = varargin{i+1};
                            i = i+2;
                        else
                            error('Invalid ''Margin'' value.');
                        end
                    case {'Start'}
                        if ischar(varargin{i+1})
                            startMethod = varargin{i+1};
                            i = i+2;
                        else
                            error('Invalid ''Start'' method.');
                        end
                    case {'Options'}
                        opts  = varargin{i+1};
                        maxIt = opts.MaxIter;
                        thres = opts.TolFun;
                        i = i+2;
                    otherwise
                        error('Unknown ''varargin{i}'' option.');
                end
            end
            
            [N, D] = size(X);

            % Initialize mixture parameters
            switch startMethod
                case {'random'}
                    P = ones(1,K)/K;
                    T = rand(K,D) * (1-2*margin) + margin; % Component means
                case {'plus'}
                    [I, T] = kmeans(X, K, 'MaxIter', 1);
                    P = hist(I, K);
                otherwise
                    error('Unsupported ''startMethod'' start method');
            end
            
            it        = 0;
            converged = false;
            postProb  = zeros(N,K); % Conditional probabilities
            firstCancelWarning = true;
            
            while ~converged && it < maxIt
                % Compute posterior sample probabilities
                for k = 1:K
                    postProb(:,k) = prod(bsxfun(@times, T(k,:), X) ...
                         + bsxfun(@times, 1-T(k,:), ~X), 2);
                end

                % Expectation phase
                Z = bsxfun(@times, postProb, P);
                % Fix cancellation issues
                cancelled = sum(Z,2) < 10^-300;
                if any(cancelled)
                    if firstCancelWarning
                        warning('Skipping samples with probability 0.');
                        firstCancelWarning = false;
                    end
                    Z = Z(~cancelled, :);
                end
                Z = bsxfun(@rdivide, Z, sum(Z,2));

                % Maximization phase
                P = mean(Z);
                for k = 1:K
                     T(k,:) = sum(bsxfun(@times, X(~cancelled,:), ... 
                         Z(:,k)), 1) / sum(Z(:,k));
                end

                % Stop or continue iterations
                if mod(it, 10) == 0
                    obj = bmdistribution(T, P);
                    if 1 - obj.pdf(X) < thres
                        return
                    end
                end
                it = it + 1;
            end

            warning('Returned before convergence.')
            obj = bmdistribution(T, P);
            obj.NumIterations = it;
        end % function fit(X, K, varargin)
        
    end
end