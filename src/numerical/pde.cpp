#include "pde.h"
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

DiffusionEquation::DiffusionEquation(bool isDisplay_,
                                        int periodDisplay_,
                                        double tMax_,
                                        std::size_t nTime_,
                                        std::vector< std::size_t >& nMesh_,
                                        std::vector< double >& lengths_,
                                        std::function< double (std::vector< double >) >& initialCondition_,
                                        std::function< double (std::vector< double >) >& dirichlet_,
                                        std::function< double (std::vector< double >) >& neumann_,
                                        std::function< double (std::vector< double >) >& thermalCond_,
                                        std::vector< std::function< double (std::vector< double >)> >& derivThermalCond_,
                                        std::vector< std::function< double (std::vector< double >)> >& boundValue_,
                                        double lambda_,
                                        double tol_,
                                        std::size_t maxIter_
                                        )
    : isDisplay( isDisplay_), periodDisplay( periodDisplay_ ), tMax( tMax_), nTime( nTime_), dim ( nMesh_.size() ), nMesh( nMesh_ ), lengths( lengths_ ), initialCondition( initialCondition_ ), dirichlet( dirichlet_ ), neumann( neumann_ ), thermalCond( thermalCond_ ), derivThermalCond( derivThermalCond_ ), boundValue( boundValue_ ), results( nMesh ), lambda( lambda_ ), lambdaMinus( 1.0 - lambda ), tol( tol_ ), maxIter( maxIter_ )
{
    totalSize = results.getTotalSize();
    dt = tMax / double( nTime );

    dMesh.resize( dim );
    positions.resize( dim );
    for(int i = 0; i < dim; i++)
    {
        dMesh.at(i) = lengths.at(i) / double( nMesh.at(i) - 3 );
        positions.at(i).resize( nMesh.at(i) );
        for(int j = 0; j < nMesh.at(i); j++)
        {
            positions.at(i).at(j) = dMesh.at(i) * double( j - 1 );
        }
    }

    for ( std::size_t offset = 0; offset < results.getTotalSize(); offset++ )
    {
        std::vector< std::size_t > indices = results.getIndices( offset );
        std::vector< double > xHere = indicesToX( indices );
        results[offset] = initialCondition( xHere );
    }
}

MultiDimMesh<double>& DiffusionEquation::calcExplicitMethod()
{
    for (tempiT = 0; tempiT < nTime; tempiT++)
    {
        if ( isDisplay && tempiT % periodDisplay == 0) { displayResults(); }
        calcExplicitMethodStep();
    }
    return results;
}

MultiDimMesh<double>& DiffusionEquation::calcExplicitMethodStep()
{
    MultiDimMesh<double> oldResults = results;
    for ( std::size_t offset = 0; offset < totalSize; offset++ )
    {
        std::vector< std::size_t > indices = results.getIndices( offset );
        std::vector< double > xHere = indicesToX( indices );
        int minOuterIndex = results.containsIndex( indices , 0 );
        int maxOuterIndex = results.containsMaxMinusIndex( indices , 0 );
        int minEdgeIndex = results.containsIndex( indices , 1 );
        int maxEdgeIndex = results.containsMaxMinusIndex( indices , 1 );
            
        if ( minOuterIndex != -1 || maxOuterIndex != -1 )
        {
            int outerIndex;
            double factor;
            std::vector< std::size_t > indicesOne( indices.size() );
            std::vector< std::size_t > indicesTwo( indices.size() );
            std::copy(indices.begin(), indices.end(), indicesOne.begin());
            std::copy(indices.begin(), indices.end(), indicesTwo.begin());
            if ( minOuterIndex != -1 )
            {
                factor = 1.0;
                outerIndex = minOuterIndex;
                indicesOne.at(outerIndex) += 1;
                indicesTwo.at(outerIndex) += 2;
            }
            else
            {
                factor = -1.0;
                outerIndex = maxOuterIndex;
                indicesOne.at(outerIndex) -= 1;
                indicesTwo.at(outerIndex) -= 2;
            }
            std::vector< double > xEdge = indicesToX( indicesOne );
            double neumannHere = neumann( xEdge );
            double dirichletHere = dirichlet( xEdge );
            double boundHere = boundValue.at(outerIndex)( xEdge );
            if( neumannHere != 0.0 )
            {
                results[indices] = (factor * dMesh.at(outerIndex) / neumannHere) * (dirichletHere * oldResults[indicesOne] - boundHere) + oldResults[indicesTwo];
            }
            continue;
        }
        if ( minEdgeIndex != -1 || maxEdgeIndex != -1 )
        {
            int edgeIndex = std::max(minEdgeIndex, maxEdgeIndex);
            double neumannHere = neumann( xHere );
            double dirichletHere = dirichlet( xHere );
            double boundHere = boundValue.at(edgeIndex)( xHere );
            if ( neumannHere == 0.0 )
            {
                results[indices] = boundHere / dirichletHere;
                continue;
            }
        }

        for ( std::size_t iX = 0; iX < dim; iX++ )
        {
            std::vector< std::size_t > indicesPlus( indices.size() );
            std::vector< std::size_t > indicesMinus( indices.size() );
            std::copy(indices.begin(), indices.end(), indicesPlus.begin());
            std::copy(indices.begin(), indices.end(), indicesMinus.begin());
            indicesPlus.at(iX) += 1;
            indicesMinus.at(iX) -= 1;
            double oldPlus = oldResults[indicesPlus];
            double oldHere = oldResults[indices];
            double oldMinus= oldResults[indicesMinus];
            results[indices] += (dt / dMesh.at(iX)) * (derivThermalCond.at(iX)( xHere ) * ( oldPlus - oldMinus ) + thermalCond( xHere ) * ( oldPlus - 2.0 * oldHere + oldMinus) / dMesh.at(iX));
        }
    }
    return results;
}

MultiDimMesh<double>& DiffusionEquation::calcCrankNicolsonMethod()
{
    for (tempiT = 0; tempiT < nTime; tempiT++)
    {
        if ( isDisplay && tempiT % periodDisplay == 0) { displayResults(); }
        calcCrankNicolsonMethodStep();
    }
    return results;
}

MultiDimMesh<double>& DiffusionEquation::calcCrankNicolsonMethodStep()
{
    Eigen::SparseMatrix<double> matSol(totalSize, totalSize);
    Eigen::VectorXd vecSol = Eigen::VectorXd::Zero(totalSize);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.setTolerance(tol);
    solver.setMaxIterations(maxIter);

    for ( std::size_t offset = 0; offset < totalSize; offset++ )
    {
        std::vector< std::size_t > indices = results.getIndices( offset );
        std::vector< double > xHere = indicesToX( indices );
        int minOuterIndex = results.containsIndex( indices , 0 );
        int maxOuterIndex = results.containsMaxMinusIndex( indices , 0 );
        int minEdgeIndex = results.containsIndex( indices , 1 );
        int maxEdgeIndex = results.containsMaxMinusIndex( indices , 1 );
            
        if ( minOuterIndex != -1 || maxOuterIndex != -1 )
        {
            int outerIndex;
            double factor;
            std::vector< std::size_t > indicesOne( indices );
            std::vector< std::size_t > indicesTwo( indices );
            if ( minOuterIndex != -1 )
            {
                factor = 1.0;
                outerIndex = minOuterIndex;
                indicesOne.at(outerIndex) += 1;
                indicesTwo.at(outerIndex) += 2;
            }
            else
            {
                factor = -1.0;
                outerIndex = maxOuterIndex;
                indicesOne.at(outerIndex) -= 1;
                indicesTwo.at(outerIndex) -= 2;
            }
            std::vector< double > xEdge = indicesToX( indicesOne );
            double neumannHere = neumann( xEdge );
            if( neumannHere != 0.0 )
            {
                double dirichletHere = dirichlet( xEdge );
                std::size_t offsetOne = results.getOffset(indicesOne);
                std::size_t offsetTwo = results.getOffset(indicesTwo);
                matSol.insert(offset, offset) = factor * neumannHere / dMesh.at(outerIndex);
                matSol.insert(offset, offsetOne) = dirichletHere;
                matSol.insert(offset, offsetTwo) = - factor * neumannHere / dMesh.at(outerIndex);
                vecSol(offset) = boundValue.at(outerIndex)( xEdge );
            }
            else
            {
                matSol.insert(offset, offset) = 1.0;
                matSol.insert(offset, results.getOffset(indicesOne)) = -1.0;
            }
            continue;
        }
        if ( minEdgeIndex != -1 || maxEdgeIndex != -1 )
        {
            int edgeIndex = std::max(minEdgeIndex, maxEdgeIndex);
            if ( neumann( xHere ) == 0.0 )
            {
                matSol.insert(offset, offset) = dirichlet( xHere );
                vecSol(offset) = boundValue.at(edgeIndex)( xHere );
                continue;
            }
        }

        double matElem = 1.0;
        double oldHere = results[offset];
        vecSol(offset) = oldHere;
        for ( std::size_t iX = 0; iX < dim; iX++ )
        {
            std::vector< std::size_t > indicesPlus( indices );
            std::vector< std::size_t > indicesMinus( indices );
            indicesPlus.at(iX) += 1;
            indicesMinus.at(iX) -= 1;

            double temp1 = derivThermalCond.at(iX)( xHere ) * dt / ( 2.0 * dMesh.at(iX) );
            double temp2 = thermalCond( xHere ) * dt / ( dMesh.at(iX) * dMesh.at(iX) );

            std::size_t offsetPlus = results.getOffset(indicesPlus);
            std::size_t offsetMinus = results.getOffset(indicesMinus);

            double oldPlus = results[offsetPlus];
            double oldMinus= results[offsetMinus];

            matElem += 2.0 * temp2 * lambda;
            matSol.insert(offset, offsetPlus) = - ( temp1 + temp2 ) * lambda;
            matSol.insert(offset, offsetMinus) = ( temp1 - temp2 ) * lambda;
            vecSol(offset) += (temp1 * (oldPlus - oldMinus) + temp2 * (oldPlus - 2.0 * oldHere + oldMinus)) * lambdaMinus;
        }
        matSol.insert(offset, offset) = matElem;
    }
    
    solver.compute(matSol);
    Eigen::VectorXd resultVector = solver.solve(vecSol);
    

    for ( std::size_t offset = 0; offset < totalSize; offset++ )
    {
        results[offset] = resultVector(offset);
    }
    return results;
}

void DiffusionEquation::displayResults()
{
    double tempT = tempiT * dt;
    for ( std::size_t offset = 0; offset < results.getTotalSize(); offset++ )
    {
        std::cout << tempT << ",";
        std::vector< std::size_t > indices = results.getIndices( offset );
        std::vector< double > xHere = indicesToX( indices );
        for ( std::size_t iX = 0; iX < dim; iX++ )
        {
            std::cout << xHere.at(iX) << ",";
        }
        std::cout << results[offset] << std::endl;
    }
}
