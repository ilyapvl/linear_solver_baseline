#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include <string>


double f1(double x) { return std::sin(x * x); }
double df1(double x) { return 2.0 * x * std::cos(x * x); }

double f2(double x) { return std::cos(std::sin(x)); }
double df2(double x) { return -std::sin(std::sin(x)) * std::cos(x); }

double f3(double x) { return std::exp(std::sin(std::cos(x))); }
double df3(double x) { return -std::sin(x) * std::cos(std::cos(x)) * std::exp(std::sin(std::cos(x))); }

double f4(double x) { return std::log(x + 3.0); }
double df4(double x) { return 1.0 / (x + 3.0); }

double f5(double x) { return std::sqrt(x + 3.0); }
double df5(double x) { return 1.0 / (2.0 * std::sqrt(x + 3.0)); }


double method1(const std::function<double(double)>& f, double x, double h)
{
    return (f(x + h) - f(x)) / h;
}

double method2(const std::function<double(double)>& f, double x, double h)
{
    return (f(x) - f(x - h)) / h;
}

double method3(const std::function<double(double)>& f, double x, double h)
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double method4(const std::function<double(double)>& f, double x, double h)
{
    return (4.0 / 3.0) * (f(x + h) - f(x - h)) / (2.0 * h) 
         - (1.0 / 3.0) * (f(x + 2.0 * h) - f(x - 2.0 * h)) / (4.0 * h);
}

double method5(const std::function<double(double)>& f, double x, double h)
{
    return (3.0 / 2.0) * (f(x + h) - f(x - h)) / (2.0 * h) 
         - (3.0 / 5.0) * (f(x + 2.0 * h) - f(x - 2.0 * h)) / (4.0 * h) 
         + (1.0 / 10.0) * (f(x + 3.0 * h) - f(x - 3.0 * h)) / (6.0 * h);
}

int main(const int argc, char* argv[])
{

    const double x0 = std::atoi(argv[1]);

    std::vector<std::pair<std::string, std::function<double(double)>>> funcs = {
        {"func1_sin_x2", f1},
        {"func2_cos_sin_x", f2},
        {"func3_exp_sin_cos_x", f3},
        {"func4_ln_x_plus_3", f4},
        {"func5_sqrt_x_plus_3", f5}
    };

    std::vector<std::function<double(double)>> derivs = {df1, df2, df3, df4, df5};
    
    std::vector<std::function<double(const std::function<double(double)>&, double, double)>> methods = {
        method1, method2, method3, method4, method5
    };


    for (size_t i = 0; i < funcs.size(); ++i)
    {
        std::string filename = funcs[i].first + "_errors.csv";
        std::ofstream out(filename);
        





        out << "h,method1,method2,method3,method4,method5\n";

        double exact_deriv = derivs[i](x0);



        for (int n = 1; n <= 45; n++)
        {
            double h = 2.0 / std::pow(2.0, n);
            out << h;

            for (const auto& method : methods)
            {
                double approx_deriv = method(funcs[i].second, x0, h);
                double error = std::abs(approx_deriv - exact_deriv);
                out << "," << error;
            }
            out << "\n";
        }
        
        out.close();
    }

    return 0;
}
