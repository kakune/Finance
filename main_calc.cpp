#include <iostream>
#include <vector>
#include <cstddef>
#include <limits>
#include <variant>
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
        Obj->set_strike(strike);
   }
}

template <class T, class U>
void setOption(const std::shared_ptr<T>& Obj, U& SpotObj)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of Option");
    Obj->set_spot_paths(SpotObj.get_spots());
    Obj->calc_payoff_mean();
}

template <class T>
void generateSpots(T& SpotObj)
{
    SpotObj.generate_spots();
}

template <class T>
void setSpotInitial(T& SpotObj, Parameters& params)
{
    SpotObj.set_num_process(params.data["COMMON"]["num_process"]);
    SpotObj.set_time(params.data["COMMON"]["num_time"], params.data["COMMON"]["maturity"]);
}

void setSpot(BlackScholes& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_vol_rate(params.data[section]["vol"], params.data["COMMON"]["rate"]);
}

void setSpot(SABR& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_initial_vol(params.data[section]["initial_vol"]);
    SpotObj.set_rho_rate(params.data[section]["rho"], params.data["COMMON"]["rate"]);
    SpotObj.set_volvol_beta(params.data[section]["volvol"], params.data[section]["beta"]);
}

void setSpot(SV& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_rho_rate(params.data[section]["rho"], params.data["COMMON"]["rate"]);
    SpotObj.set_params_spot(params.data[section]["lambda"], params.data[section]["b"], params.data[section]["L"]);
    SpotObj.set_params_vol(params.data[section]["theta"], params.data[section]["eta"]);
}

void showIV(const std::vector<double>& variable, const std::vector<double>& result, Parameters& params, const std::string& section)
{
    params.initVariable();
    ImpliedVolatility ImpVolObj;
    ImpVolObj.setRate(params.data["COMMON"]["rate"]);
    ImpVolObj.setMaturity(params.data["COMMON"]["maturity"]);
    ImpVolObj.setMaxIter(params.data["COMMON"]["maxIter"]);
    ImpVolObj.setTol(params.data["COMMON"]["tol"]);
    ImpVolObj.setReg(params.data["COMMON"]["volMin"], params.data["COMMON"]["volMax"]);

    int i = 0;
    while( params.updateVariable() )
    {
        ImpVolObj.setSpot( params.data[section]["initial_spot"] );
        ImpVolObj.setStrike( params.data[section]["strike"] );
        ImpVolObj.setCallValue( result.at(i) );
        std::cout << variable.at(i) << "," << ImpVolObj.calcImpliedVol() << std::endl;
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

    std::variant<BlackScholes, SABR, SV> SpotObj;
    
    if(params.sdata[section]["Model"] == "BS" || params.sdata[section]["Model"] == "BlackScholes")
    {SpotObj = BlackScholes();}
    else if (params.sdata[section]["Model"] == "SABR")
    {SpotObj = SABR();}
    else if (params.sdata[section]["Model"] == "SV")
    {SpotObj = SV();}
    else
    {
        std::cerr << "Error : Such model has not been implemented. " << std::endl;
        std::exit(1);
    }

    std::visit([&](auto&& arg){setSpotInitial(arg, params);}, SpotObj);
    std::visit([&](auto&& arg){setSpot(arg, params, section);}, SpotObj);
    if(params.sdata["COMMON"]["path_each"] != "true"){std::visit([&](auto&& arg){generateSpots(arg);}, SpotObj);}
    while( params.updateVariable() )
    {
        setStrikeToOptions(CallVector, params.data[section]["strike"]);
        setStrikeToOptions(PutVector, params.data[section]["strike"]);

        std::visit([&](auto&& arg){setSpot(arg, params, section);}, SpotObj);
        if(params.sdata["COMMON"]["path_each"] == "true"){std::visit([&](auto&& arg){generateSpots(arg);}, SpotObj);}
        std::visit([&](auto&& arg){setOption(OptionObj, arg);}, SpotObj);

        variable.push_back(params.getVariable());
        result.push_back(OptionObj->get_payoff_mean() * df);
    }
    }
    else if (params.sdata[section]["method"] == "Fourier")
    {
    FourierPricingSV FSVObj;
    while( params.updateVariable() )
    {
        FSVObj.setMaturity(params.data["COMMON"]["maturity"]);
        FSVObj.setRate(params.data["COMMON"]["rate"]);
        FSVObj.setRho(params.data[section]["rho"]);
        FSVObj.setParamsSpot(params.data[section]["lambda"],params.data[section]["b"],params.data[section]["L"]);
        FSVObj.setParamsVol(params.data[section]["theta"],params.data[section]["eta"]);
        FSVObj.setSpot(params.data[section]["initial_spot"]);
        FSVObj.setStrike(params.data[section]["strike"]);

        variable.push_back(params.getVariable());
        result.push_back(FSVObj.calcCallValue() * df);
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
