#ifndef NUMERICAL_PDE_H
#define NUMERICAL_PDE_H

#include <vector>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <iostream>
template <typename T>
class MultiDimMesh
{
private:
    std::size_t dim, totalSize;
    std::vector< std::size_t > nMesh;
    std::vector< T > data;

public:
    MultiDimMesh(const std::vector< std::size_t >& nMesh_) 
        : nMesh( nMesh_ ), dim( nMesh_.size() )
    {
        totalSize = std::accumulate(nMesh.begin(), nMesh.end(), 
                            1, std::multiplies<size_t>());
        data.resize(totalSize);
    }

    size_t getOffset(const std::vector<size_t>& indices) const {
        size_t offset = 0;
        size_t stride = 1;

        if ( indices.size() != nMesh.size() )
        {
            throw std::invalid_argument("Dimension mismatch");
        }

        for ( int i = indices.size() - 1; i >= 0; --i )
        {
            if(indices[i] >= nMesh[i]) {
                throw std::out_of_range("Index out of range");
            }
            offset += indices[i] * stride;
            stride *= nMesh[i];
        }

        return offset;
    }

    std::vector< std::size_t > getIndices( std::size_t offset ) const {
        std::vector< std::size_t > indices(nMesh.size());
        for (int i = nMesh.size() - 1; i >= 0; --i) {
            indices[i] = offset % nMesh[i];
            offset /= nMesh[i];
        }
        return indices;
    }

    int containsIndex( const std::vector< std::size_t >& indices , const std::size_t& value ) const 
    {
        for ( int i = 0; i < indices.size(); ++i )
        {
            if ( indices[i] == value )
                return i;
        }
        return -1;
    }
    int containsMaxMinusIndex( const std::vector< std::size_t >& indices , const std::size_t& value ) const 
    {
        for ( int i = 0; i < indices.size(); ++i )
        {
            if ( indices[i] == nMesh[i] - 1 - value )
                return i;
        }
        return -1;
    }
    int containsIndex( const std::size_t& offset , const std::size_t& value ) const 
    {
        std::vector< std::size_t > indices = getIndices( offset );
        return containsIndex( indices, value );
    }
    int containsMaxMinusIndex( const std::size_t& offset , const std::size_t& value ) const 
    {
        std::vector< std::size_t > indices = getIndices( offset );
        return containsMaxMinusIndex( indices, value );
    }

    T& operator[](const std::vector< std::size_t >& indices)
    {
        return data[getOffset(indices)];
    }

    const T& operator[](const std::vector< std::size_t >& indices) const
    {
        return data[getOffset(indices)];
    }

    T& operator[](const std::size_t& offset)
    {
        return data[offset];
    }

    const T& operator[](const std::size_t& offset) const
    {
        return data[offset];
    }

    std::size_t getTotalSize() const
    {
        return totalSize;
    }
    std::size_t getDim() const
    {
        return dim;
    }

};
// region : [regMin.at(0), regMax.at(0)] * [regMin.at(1), regMax.at(1)] * ... 
// dy/dt = dR/dxi dy/dxi + R d^2 y / dx^2
// y(x, 0) = initialCondition( x )
// dirichlet(x) * y(x, t) + neumann(x) * dy/dxi (x, t) = bound_i(x)

class DiffusionEquation
{
private:
    bool isDisplay;
    int periodDisplay;
    std::size_t tempiT;
    double tMax, dt;
    std::size_t nTime;
    std::size_t dim;
    std::vector< std::size_t > nMesh; // To consider Neumann b.c., nMesh must include outer mesh 
    std::vector< double > dMesh, lengths;
    std::function< double (std::vector< double >) > initialCondition, dirichlet, neumann, thermalCond;
    std::vector< std::function< double (std::vector< double >)> > derivThermalCond, boundValue;
    MultiDimMesh<double> results;
    std::vector< std::vector< double > > positions;
    std::size_t totalSize;
    double lambda, lambdaMinus, tol;
    std::size_t maxIter;
public:
    DiffusionEquation(bool isDisplay_,
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
                    );
    std::vector< std::vector< double > >& getPositions() { return positions; }
    MultiDimMesh<double>& calcExplicitMethodStep();
    MultiDimMesh<double>& calcExplicitMethod();
    MultiDimMesh<double>& calcCrankNicolsonMethodStep();
    MultiDimMesh<double>& calcCrankNicolsonMethod();
    std::vector< double > indicesToX( std::vector< std::size_t > indices )
    {
        std::vector< double > res( dim, 0.0 );
        for ( std::size_t iX = 0; iX < dim; iX++ )
        {
            res.at(iX) = positions.at(iX).at(indices.at(iX));
        }
        return res;
    }
    void displayResults();
    // std::vector< std::vector< double > > getPositions {}
};

class DiffusionEquationBuilder
{
private:
    bool isDisplay = false;
    int periodDisplay = 20;
    double tMax;
    double lambda = 0.5;
    double tol = 1e-8;
    std::size_t maxIter = 100;
    std::size_t nTime;
    std::vector< std::size_t > nMesh; // To consider Neumann b.c., nMesh must include outer mesh 
    std::vector< double > lengths;
    std::function< double (std::vector< double >) > initialCondition, dirichlet, neumann, thermalCond;
    std::vector< std::function< double (std::vector< double >)> > derivThermalCond, boundValue;
public:
    DiffusionEquationBuilder& setIsDisplay( bool isDisplay_ )
    {
        isDisplay = isDisplay_;
        return *this;
    }
    DiffusionEquationBuilder& setPeriodDisplay( int periodDisplay_ )
    {
        periodDisplay = periodDisplay_;
        return *this;
    }
    DiffusionEquationBuilder& setTMax( double tMax_ )
    {
        tMax = tMax_;
        return *this;
    }
    DiffusionEquationBuilder& setNTime( std::size_t nTime_ )
    {
        nTime = nTime_;
        return *this;
    }
    DiffusionEquationBuilder& setNMesh( const std::vector< std::size_t >& nMesh_ )
    {
        nMesh = nMesh_;
        return *this;
    }
    DiffusionEquationBuilder& setNMesh( std::size_t nMesh_ )
    {
        nMesh = std::vector< std::size_t > (1, nMesh_);
        return *this;
    }
    DiffusionEquationBuilder& setLengths( const std::vector< double >& lengths_ )
    {
        lengths = lengths_;
        return *this;
    }
    DiffusionEquationBuilder& setLengths( double lengths_ )
    {
        lengths = std::vector<double>(1, lengths_);
        return *this;
    }
    DiffusionEquationBuilder& setInitialCondition( const std::function< double (std::vector< double >) >& initialCondition_ )
    {
        initialCondition = initialCondition_;
        return *this;
    }
    DiffusionEquationBuilder& setInitialCondition( const std::function< double (double) >& initialCondition_ )
    {
        initialCondition = [initialCondition_](std::vector<double> x){ return initialCondition_(x.at(0));} ;
        return *this;
    }
    DiffusionEquationBuilder& setDirichlet( const std::function< double (std::vector< double >) >& dirichlet_ )
    {
        dirichlet = dirichlet_;
        return *this;
    }
    DiffusionEquationBuilder& setDirichlet( const std::function< double (double) >& dirichlet_ )
    {
        dirichlet = [dirichlet_](std::vector<double> x){ return dirichlet_(x.at(0));} ;
        return *this;
    }
    DiffusionEquationBuilder& setNeumann( const std::function< double (std::vector< double >) >& neumann_ )
    {
        neumann = neumann_;
        return *this;
    }
    DiffusionEquationBuilder& setNeumann( const std::function< double (double) >& neumann_ )
    {
        neumann = [neumann_](std::vector<double> x){ return neumann_(x.at(0));} ;
        return *this;
    }
    DiffusionEquationBuilder& setThermalCond( const std::function< double (std::vector< double >) >& thermalCond_ )
    {
        thermalCond = thermalCond_;
        return *this;
    }
    DiffusionEquationBuilder& setThermalCond( const std::function< double (double) >& thermalCond_ )
    {
        thermalCond = [thermalCond_](std::vector<double> x){ return thermalCond_(x.at(0));};
        return *this;
    }
    DiffusionEquationBuilder& setDerivThermalCond( const std::vector< std::function< double (std::vector< double >)> >& derivThermalCond_ )
    {
        derivThermalCond = derivThermalCond_;
        return *this;
    }
    DiffusionEquationBuilder& setDerivThermalCond( const std::function< double (double) >& derivThermalCond_ )
    {
        derivThermalCond.resize(1, [derivThermalCond_](std::vector<double> x){ return derivThermalCond_(x.at(0));});
        return *this;
    }
    DiffusionEquationBuilder& setDerivThermalCond( const std::vector< std::function< double (double) > >& derivThermalCond_ )
    {
        derivThermalCond.resize(1, [derivThermalCond_](std::vector<double> x){ return derivThermalCond_.at(0)(x.at(0));} );
        return *this;
    }
    DiffusionEquationBuilder& setDerivThermalCond( const std::function< double (std::vector< double >) >& derivThermalCond_ )
    {
        derivThermalCond.resize(1, derivThermalCond_);
        return *this;
    }
    DiffusionEquationBuilder& setBoundValue( const std::vector< std::function< double (std::vector< double >)> >& boundValue_ )
    {
        boundValue = boundValue_;
        return *this;
    }
    DiffusionEquationBuilder& setBoundValue( const std::function< double (double) >& boundValue_ )
    {
        boundValue.resize(1, [boundValue_](std::vector<double> x){ return boundValue_(x.at(0));});
        return *this;
    }
    DiffusionEquationBuilder& setBoundValue( const std::vector< std::function< double (double) > >& boundValue_ )
    {
        boundValue.resize(1, [boundValue_](std::vector<double> x){ return boundValue_.at(0)(x.at(0));} );
        return *this;
    }
    DiffusionEquationBuilder& setBoundValue( const std::function< double (std::vector< double >) >& boundValue_ )
    {
        boundValue.resize(1, boundValue_);
        return *this;
    }
    DiffusionEquationBuilder& setLambda( double lambda_ )
    {
        lambda = lambda_;
        return *this;
    }
    DiffusionEquationBuilder& setTol( double tol_ )
    {
        tol = tol_;
        return *this;
    }
    DiffusionEquationBuilder& setMaxIter( std::size_t maxIter_ )
    {
        maxIter = maxIter_;
        return *this;
    }
    DiffusionEquation build()
    {
        return DiffusionEquation(isDisplay, periodDisplay, tMax, nTime, nMesh, lengths, initialCondition, dirichlet, neumann, thermalCond, derivThermalCond, boundValue, lambda, tol, maxIter);
    }

};

#endif