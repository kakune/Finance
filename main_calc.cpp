#include <iostream>
#include <vector>
#include <cstddef>
#include <limits>
#include <variant>
#include <utility>
#include <core/params.h>
#include <monte_carlo/process.h>
#include <option/european.h>
#include <option/option_core.h>
#include <option/implied_vol.h>
#include <fourier/vanilla.h>


template <class T>
void setStrikeToOptions(const std::vector<std::shared_ptr<T> >& Objs, const double& strike)
{
   for(std::shared_ptr<T> Obj : Objs)
   {
        Obj->setStrike(strike);
   }
}

template <class T>
void setOption(const std::shared_ptr<T>& Obj, const std::vector < std::vector<double> >& spots)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of Option");
    Obj->setSpotPaths(spots);
    Obj->calcPayoffMean();
}


template <class T>
void setSpotInitial(T& SpotBuilder, Parameters& params)
{
    SpotBuilder.setNumProcess( params.data["COMMON"]["num_process"]);
    SpotBuilder.setNumTime( params.data["COMMON"]["num_time"] );
    SpotBuilder.setMaturity( params.data["COMMON"]["maturity"] );
}

void setSpot(BlackScholesBuilder& SpotBuilder, Parameters& params, const std::string& section)
{
    SpotBuilder.setInitialSpot(params.data[section]["initial_spot"]);
    SpotBuilder.setVol(params.data[section]["vol"]);
    SpotBuilder.setRate(params.data["COMMON"]["rate"]);
}

void setSpot(SABRBuilder& SpotBuilder, Parameters& params, const std::string& section)
{
    SpotBuilder.setInitialSpot(params.data[section]["initial_spot"]);
    SpotBuilder.setInitialVol(params.data[section]["initial_vol"]);
    SpotBuilder.setRho(params.data[section]["rho"]);
    SpotBuilder.setRate(params.data["COMMON"]["rate"]);
    SpotBuilder.setVolvol(params.data[section]["volvol"]);
    SpotBuilder.setBeta(params.data[section]["beta"]);
}

void setSpot(SVBuilder& SpotBuilder, Parameters& params, const std::string& section)
{
    SpotBuilder.setInitialSpot(params.data[section]["initial_spot"]);
    SpotBuilder.setRho(params.data[section]["rho"]);
    SpotBuilder.setRate(params.data["COMMON"]["rate"]);
    SpotBuilder.setLambda(params.data[section]["lambda"]);
    SpotBuilder.setB(params.data[section]["b"]);
    SpotBuilder.setL(params.data[section]["L"]);
    SpotBuilder.setTheta(params.data[section]["theta"]);
    SpotBuilder.setEta(params.data[section]["eta"]);
}

void showIV(const std::vector<double>& variable, const std::vector<double>& result, Parameters& params, const std::string& section)
{
    params.initVariable();
    ImpliedVolatilityBuilder ImpVolBuilder;
    ImpVolBuilder.setRate(params.data["COMMON"]["rate"]);
    ImpVolBuilder.setMaturity(params.data["COMMON"]["maturity"]);
    ImpVolBuilder.setMaxIter(params.data["COMMON"]["maxIter"]);
    ImpVolBuilder.setTol(params.data["COMMON"]["tol"]);
    ImpVolBuilder.setVolReg(params.data["COMMON"]["volMin"], params.data["COMMON"]["volMax"]);

    int i = 0;
    while( params.updateVariable() )
    {
        ImpVolBuilder.setSpot( params.data[section]["initial_spot"] );
        ImpVolBuilder.setStrike( params.data[section]["strike"] );
        std::cout << variable.at(i) << "," << ImpVolBuilder.build().calcImpliedVol( result.at(i) ) << std::endl;
        i++;
    }
}

void showResult(const std::vector<double>& variable, const std::vector<double>& result)
{
    for(int i = 0; i < variable.size(); i++)
    {
        std::cout << variable.at(i) << "," << result.at(i) << std::endl;
    }
}


int main(int argc, char* argv[])
{
    Parameters params;
    params.readParameters(argv[1]);
    std::string section = argv[2];

    params.setVariable("VARIABLE", "variable", section);
    std::vector<double> variable, result;
    double df = std::exp(- params.data["COMMON"]["maturity"] * params.data["COMMON"]["rate"]);

    if (params.sdata[section]["method"] == "MC")
    {
        std::shared_ptr<EuropeanCallOption> OptionObj = std::make_shared<EuropeanCallOption>();
        std::shared_ptr<EuropeanPutOption> PutObj = std::make_shared<EuropeanPutOption>();
        std::vector< std::shared_ptr<Call> > CallVector;
        std::vector< std::shared_ptr<Put> > PutVector;
        CallVector.push_back(OptionObj);
        PutVector.push_back(PutObj);
        // PutVector.push_back(PutObj);

        // std::shared_ptr<VirtualOption> OptionObj = CallObj + PutObj;

        std::variant<BlackScholesBuilder, SABRBuilder, SVBuilder> SpotBuilder;
        std::vector< std::vector<double> > spots;

        if(params.sdata[section]["Model"] == "BS" || params.sdata[section]["Model"] == "BlackScholes")
        {SpotBuilder = BlackScholesBuilder();}
        else if (params.sdata[section]["Model"] == "SABR")
        {SpotBuilder = SABRBuilder();}
        else if (params.sdata[section]["Model"] == "SV")
        {SpotBuilder = SVBuilder();}
        else
        {
            std::cerr << "Error : Such model has not been implemented. " << std::endl;
            std::exit(1);
        }

        std::visit([&](auto&& arg){setSpotInitial(arg, params);}, SpotBuilder);
        std::visit([&](auto&& arg){setSpot(arg, params, section);}, SpotBuilder);
        if ( params.sdata["COMMON"]["path_each"] != "true" )
        {
            std::visit([&](auto&& builder)
            {
                spots = std::move(builder.build().generateSpots());
            }, SpotBuilder);
        }
        while( params.updateVariable() )
        {
            setStrikeToOptions(CallVector, params.data[section]["strike"]);
            setStrikeToOptions(PutVector, params.data[section]["strike"]);

            std::visit([&](auto&& arg){setSpot(arg, params, section);}, SpotBuilder);
            if(params.sdata["COMMON"]["path_each"] == "true")
            {
                std::visit([&](auto&& builder)
                {
                    spots = std::move(builder.build().generateSpots());
                }, SpotBuilder);
            }
            setOption(OptionObj, spots);

            variable.push_back(params.getVariable());
            result.push_back(OptionObj->getPayoffMean() * df);
        }
    }
    else if (params.sdata[section]["method"] == "Fourier")
    {
        FourierPricingSVBuilder FSVBuilder;
        while( params.updateVariable() )
        {
            FSVBuilder.setMaturity( params.data["COMMON"]["maturity"] );
            FSVBuilder.setRate( params.data["COMMON"]["rate"] );
            FSVBuilder.setRho( params.data[section]["rho"] );
            FSVBuilder.setLambda( params.data[section]["lambda"] );
            FSVBuilder.setB( params.data[section]["b"] );
            FSVBuilder.setL( params.data[section]["L"] );
            FSVBuilder.setTheta( params.data[section]["theta"]);
            FSVBuilder.setEta( params.data[section]["eta"] );
            FSVBuilder.setSpot( params.data[section]["initial_spot"] );
            FSVBuilder.setStrike( params.data[section]["strike"] );

            variable.push_back(params.getVariable());
            result.push_back(FSVBuilder.build().calcCallValue() * df);
        }
    }
    if ( params.sdata["COMMON"]["objective"] == "IV" || params.sdata["COMMON"]["objective"] == "impliedVolatility" )
    {
        showIV(variable, result, params, section);
    }
    else
    {
        showResult(variable, result);
    }

    return 0;
}
