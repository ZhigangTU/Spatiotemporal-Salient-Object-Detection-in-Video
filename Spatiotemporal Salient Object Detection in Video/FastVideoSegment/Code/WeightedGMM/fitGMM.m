% Function to fit a GMM to points weighted by importance
%
%    Copyright (C) 2013  Anestis Papazoglou
%
%    You can redistribute and/or modify this software for non-commercial use
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    For commercial use, contact the author for licensing options.
%
%    Contact: a.papazoglou@sms.ed.ac.uk

function [model, MaxSglarLab] = fitGMM( clusters, data, dataWeights )

    if( ~exist( 'dataWeights', 'var' ) )
        % TODO: Classic GMM, all points are equally important
    else
        if( isinteger( data ) )
            data = double( data );
        end
            
        labels = vectorQuantize( data', clusters );
        
        [mu, sigma, weights, MaxSglarLab] = eStep(data, dataWeights, labelsToGamma( labels, clusters ), clusters);
        gamma = mStep( data, mu, sigma, weights );
        
        oldLikelihood = 10;
        likelihood = 5;
        while( abs( oldLikelihood - likelihood ) < 1e-5 )
            % Fit gaussians
            [ mu, sigma, weights ] = eStep( data, dataWeights, ...
                gamma, clusters );
            % Recompute responsibilities
            gamma = mStep( data, mu, sigma, weights );
            % Compute likelihood
            oldLikelihood = likelihood;
            likelihood = logLikelihood( data, dataWeights, mu, sigma, weights );
        end
        
        model = gmdistribution( mu, sigma, weights );
    end

end

function [mu, sigma, weights, MaxSglarLab] = eStep( data, dataWeights, gamma, clusters )

    epss = eps( 'single' );
    
    MaxSglarNum = 0;
    
    points = size( data, 1 );
    dimensions = size( data, 2 );
    
    sigma = zeros( dimensions, dimensions, clusters );
    
    dWeights = sparse( 1: points, 1: points, double( dataWeights ) );
    Gamma = gamma * dWeights;
    normaliser = ...
        sparse( 1: clusters, 1: clusters, 1 ./ sum( Gamma, 2 ) );
    
    mu = normaliser * ( ( Gamma * data ) );
    
    for cluster = 1: clusters
        dataMinusMu = bsxfun( @minus, data, mu( cluster, : ) );
        clusterSigma = normaliser( cluster, cluster ) * ( bsxfun( ...
            @times, dataMinusMu, Gamma( cluster, : )' )  )' * dataMinusMu;
        % Copy upper triangular part of sigma to lower part to ensure that
        % there are no floating precision errors that make the matrix asymmetrical
        clusterSigma( tril( true( dimensions ), -1 ) ) = ...
            clusterSigma( triu( true( dimensions ), 1 ) );
        sigma( :, :, cluster ) = clusterSigma;

        % Check that sigma is well conditioned
        if( rcond( sigma( :, :, cluster ) ) < epss )
            warning( 'fitGMM:illConditionedSigma', 'Warning, degenerate cluster detected' );
            sigma( :, :, cluster ) = sigma( :, :, cluster ) + 1e-6 * eye( dimensions );
            MaxSglarNum = MaxSglarNum+1;
        end
    end
    
    if MaxSglarNum>0.6*clusters
        MaxSglarLab = true;
    else
        MaxSglarLab = false;
    end

    weights = sum( Gamma, 2 ) / sum( dataWeights );
    
end

function gamma = mStep( data, mu, sigma, weights )

    points = size( data, 1 );
    clusters = length( weights );
    
    gamma = zeros( clusters, points );
    normaliser = zeros( 1, points );
    for cluster = 1: clusters

        gamma( cluster, : ) = weights( cluster ) * mvnpdf( data, mu( cluster, : ), sigma( :, :, cluster ) );

        normaliser = normaliser + gamma( cluster, : );
    end

    gamma = bsxfun( @rdivide, gamma, normaliser );

end

function likelihood = logLikelihood( data, dataWeights, mu, sigma, ...
    weights )

    points = size( data, 1 );
    clusters = length( weights );

    dataLikelihood = zeros( points, 1 );
    for( cluster = 1: clusters )
        dataLikelihood = dataLikelihood + weights( cluster ) * ...
            mvnpdf( data, mu( cluster, : ), sigma( :, :, cluster ) );
    end

    dataLikelihood = log( dataLikelihood );
    likelihood = sum( dataWeights .* dataLikelihood );
    
end

function gamma = labelsToGamma( labels, clusters )

    gamma = zeros( clusters, length( labels ) );
    for( i = 1: clusters )
        gamma( i, : ) = labels == i;
    end
    
end
