#include <iostream>
#include <vector>
#include <cstddef>
#include <limits>
#include <variant>
#include <core/params.h>
#include <monte_carlo/process.h>
#include <option/european.h>
#include <option/option_core.h>
#include <numerical/find_root.h>
#include <option/implied_vol.h>


template <class T>
void set_strike_to_Options(const std::vector<std::shared_ptr<T> >& Objs, const double& strike)
{
   for(std::shared_ptr<T> Obj : Objs)
   {
        Obj->set_strike(strike);
   }
}

template <class T, class U>
void set_Option(const std::shared_ptr<T>& Obj, U& SpotObj)
{
    static_assert(std::is_base_of<Option, T>::value, "T must be a descendant of Option");
    Obj->set_spot_paths(SpotObj.get_spots());
    Obj->calc_payoff_mean();
}

template <class T>
void set_Spot_initial(T& SpotObj, Parameters& params)
{
    SpotObj.set_num_process(params.data["COMMON"]["num_process"]);
    SpotObj.set_time(params.data["COMMON"]["num_time"], params.data["COMMON"]["maturity"]);
}

void set_Spot(BlackScholes& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_vol_rate(params.data[section]["vol"], params.data["COMMON"]["rate"]);
    SpotObj.generate_spots();
}

void set_Spot(SABR& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_initial_vol(params.data[section]["initial_vol"]);
    SpotObj.set_rho_rate(params.data[section]["rho"], params.data["COMMON"]["rate"]);
    SpotObj.set_volvol_beta(params.data[section]["volvol"], params.data[section]["beta"]);
    SpotObj.generate_spots();
}

void set_Spot(SV& SpotObj, Parameters& params, const std::string& section)
{
    SpotObj.set_initial_spot(params.data[section]["initial_spot"]);
    SpotObj.set_rho_rate(params.data[section]["rho"], params.data["COMMON"]["rate"]);
    SpotObj.set_params_spot(params.data[section]["lambda"], params.data[section]["b"], params.data[section]["L"]);
    SpotObj.set_params_vol(params.data[section]["theta"], params.data[section]["eta"]);
    SpotObj.generate_spots();
}


void set_BS_Analytical(AnalyticalBlackScholesCall& Obj, Parameters& params, const std::string& section)
{
    Obj.set_spot(params.data[section]["initial_spot"]);
    Obj.set_strike(params.data[section]["strike"]);
    Obj.set_vol(params.data[section]["vol"]);
    Obj.set_rate(params.data[section]["rate"]);
    Obj.set_maturity(params.data["COMMON"]["maturity"]);
}

int main(int argc, char* argv[]) {
    Parameters params;
    params.readParameters(argv[1]);
    std::string section = argv[2];

    params.setVariable("VARIABLE", "variable", section);

    std::vector<double> variable, result;
    std::shared_ptr<EuropeanCallOption> OptionObj = std::make_shared<EuropeanCallOption>();
    std::shared_ptr<EuropeanPutOption> PutObj = std::make_shared<EuropeanPutOption>();
    std::vector< std::shared_ptr<Call> > CallVector;
    std::vector< std::shared_ptr<Put> > PutVector;
    CallVector.push_back(OptionObj);
    PutVector.push_back(PutObj);

    double df = std::exp(- params.data["COMMON"]["maturity"] * params.data["COMMON"]["rate"]);
    // PutVector.push_back(PutObj);

    // std::shared_ptr<VirtualOption> OptionObj = CallObj + PutObj;

    std::variant<BlackScholes, SABR, SV> SpotObj;
    
    if(params.sdata[section]["Model"] == "BS" || params.sdata[section]["Model"] == "BlackScholes")
    {
        SpotObj = BlackScholes();
    }
    else if (params.sdata[section]["Model"] == "SABR")
    {
        SpotObj = SABR();
    }
    else if (params.sdata[section]["Model"] == "SV")
    {
        SpotObj = SV();
    }
    else
    {
        std::cerr << "Error : Such model has not been implemented. " << std::endl;
        std::exit(1);
    }

    std::visit([&](auto&& arg){set_Spot_initial(arg, params);}, SpotObj);
    if(params.sdata["COMMON"]["path_each"] != "true"){std::visit([&](auto&& arg){set_Spot(arg, params, section);}, SpotObj);}
    while(params.updateVariable())
    {
        set_strike_to_Options(CallVector, params.data[section]["strike"]);
        set_strike_to_Options(PutVector, params.data[section]["strike"]);

        if(params.sdata["COMMON"]["path_each"] == "true"){std::visit([&](auto&& arg){set_Spot(arg, params, section);}, SpotObj);}
        std::visit([&](auto&& arg){set_Option(OptionObj, arg);}, SpotObj);

        variable.push_back(params.getVariable());
        result.push_back(OptionObj->get_payoff_mean() * df);
    }



    double reg_min = 0.0001;   // Lower guess for implied volatility
    double reg_max = 100.0;    // Upper guess for implied volatility
    BrentSolver solver(reg_min, reg_max);   // Construct the Brent solver
    AnalyticalBlackScholesCall BSModel;
    set_BS_Analytical(BSModel, params, section);

    for(int i = 0; i < variable.size(); i++)
    {
        BSModel.set_strike(variable.at(i));
        // Define the function for which we want to find the root
        auto func = [&](double sigma){ BSModel.set_vol(sigma); return BSModel.payoff() - result.at(i); };
        // Find the implied volatility
        solver.set_reg(reg_min, reg_max);
        double impliedVol = solver.solve(func);
        std::cout << variable.at(i) << "," << impliedVol << std::endl;
        // std::cout << variable.at(i) << "," << result.at(i) << std::endl;
        
    }

    return 0;
}
