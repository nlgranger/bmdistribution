% Author  : Nicolas Granger <nicolas.granger@telecom-sudparis.eu>
% License : MIT (see LICENSE file)
% Source  : https://github.com/pixelou/bmdistribution

function bmm = fitbmdist(X, K, varargin)
    %FITBMDIST Fits a Binary Mixture Model to a distribution of data.
    %   bmm = FITBMDIST(X, K) Fits Binary (Bernoulli) mixture distribution to
    %   the data distribution. Rows of X correspnd to observations; columns
    %   correspond to variables. FITBMDIST fits the model by maximum 
    %   likelihood, using the Expectation-Maximization (EM) algorithm.
    %
    %   FITBMDIST treats NaNs as missing data. Rows of X with NaNs are excluded 
    %   from the fit.
    %   
    %   bmm = FITBMDIST(X, K, 'PARAM1', val1, 'PARAM2', val2, ...) allows you 
    %   to specify optional parameter name/value pairs to specify details of 
    %   the model and to control the iterative EM algorithm used to fit the 
    %   model.
    %   Parameters are:
    %   
    %       'Start'     The method used to choose initial means for the 
    %                   Bernoulli components. Choices are:
    % 
    %                   'random', to generate random means in range [.1 .9]
    %                   and uniform mixing proportion.
    %
    %                   'plus', to use the centroids of a kmeans algorithm
    %                   as the inital means. Only one iteration of kmeans
    %                   is done, and the kmeans++ initialization is used. 
    %                   Initial mixing proportion is set to each kmeans
    %                   class proportion.
    %
    %       'Options'   Options structure for the iterative EM algorithm.
    %                   The following parameters are used:
    %
    %                   'MaxIter'  Maximum number of iterations allowed.
    %                              Default is 500.
    %                   'TolFun'   Positive number giving the termination
    %                              tolerance for the log-likelyhood
    %                              function. Default is 1e-6.
    %   
    %   Reference :
    %   JUAN, Alfons et VIDAL, Enrique. Bernoulli Mixture Models for Binary 
    %   Images. In : ICPR (3). 2004. p. 367-370.
    
	bmm = bmdistribution.fit(X, K, varargin);
end