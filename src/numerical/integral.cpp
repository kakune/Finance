#include "integral.h"

double OneDimIntegrationGSL::cquadIntegral()
{
    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc(100); 
    gsl_integration_cquad(         
        &F,      //被積分関数の定義
        regMin, regMax,  //積分範囲(a,b) 
        0, 1e-10, //収束させる推定絶対誤差と相対誤差
        w,       //作業領域のサイズと 作業領域
        &result, //計算結果取得
        &error,  //推定絶対誤差取得
        &nevals  //評価回数取得
    );
    gsl_integration_cquad_workspace_free(w);
    return result;
}
double OneDimIntegrationGSL::qagiIntegral()
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100); 
    gsl_integration_qagi(         
        &F,      //被積分関数の定義
        1e-6, 1e-6, //収束させる推定絶対誤差と相対誤差
        100000,
        w,       //作業領域のサイズと 作業領域
        &result, //計算結果取得
        &error  //推定絶対誤差取得
    );
    gsl_integration_workspace_free(w);
    return result;
}
double OneDimIntegrationGSL::qagiuIntegral()
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100); 
    gsl_integration_qagiu(       
        &F,      //被積分関数の定義
        regMin,
        0, 1e-10, //収束させる推定絶対誤差と相対誤差
        10000000,
        w,       //作業領域のサイズと 作業領域
        &result, //計算結果取得
        &error  //推定絶対誤差取得
    );
    gsl_integration_workspace_free(w);
    return result;
}
double OneDimIntegrationGSL::qagilIntegral()
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100); 
    gsl_integration_qagil(         
        &F,      //被積分関数の定義
        regMax,
        1e-6, 1e-6, //収束させる推定絶対誤差と相対誤差
        100000,
        w,       //作業領域のサイズと 作業領域
        &result, //計算結果取得
        &error  //推定絶対誤差取得
    );
    gsl_integration_workspace_free(w);
    return result;
}

double oneDimIntegralCquadGSL(const std::function<double(double)>& integrand, double regMin, double regMax)
{
    OneDimIntegrationGSL GSLObj(integrand, regMin, regMax);
    return GSLObj.cquadIntegral();
}

double oneDimIntegralQagiGSL(const std::function<double(double)>& integrand)
{
    OneDimIntegrationGSL GSLObj(integrand, 0.0, 0.0);
    return GSLObj.qagiIntegral();
}

double oneDimIntegralQagiuGSL(const std::function<double(double)>& integrand, double regMin)
{
    OneDimIntegrationGSL GSLObj(integrand, regMin, 0.0);
    return GSLObj.qagiuIntegral();
}

double oneDimIntegralQagilGSL(const std::function<double(double)>& integrand, double regMax)
{
    OneDimIntegrationGSL GSLObj(integrand, 0.0, regMax);
    return GSLObj.qagilIntegral();
}

// #include <iostream>
// int main()
// {
//     std::function<double(double)> integ = [](double x){return x*x;};
//     // OneDimIntegrationGSL GSLObj;
//     // GSLObj.setReg(0.0, 1.0);
//     // GSLObj.setIntegrand(integ);
//     // GSLObj.calcIntegral();
//     std::cout << oneDimIntegral(integ, 0.0, 1.0, 50000) << std::endl;
// }