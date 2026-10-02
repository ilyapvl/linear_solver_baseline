#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <algorithm>

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

constexpr double EPS = 1e-15;

static double norm2(const Vec& v)
{
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}

static Vec matvec(const Mat& A, const Vec& x)
{
    const int n = static_cast<int>(A.size());
    Vec y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += A[i][j] * x[j];
        y[i] = s;
    }
    return y;
}

static double dot(const Vec& a, const Vec& b)
{
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
    return s;
}

static Vec sub(const Vec& a, const Vec& b)
{
    Vec r(a.size());
    for (std::size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
    return r;
}

static double residual_norm(const Mat& A, const Vec& x, const Vec& b)
{
    return norm2(sub(matvec(A, x), b));
}

static void write_history(const std::string& fname, const std::vector<double>& h)
{
    std::ofstream f(fname);
    f << std::scientific << std::setprecision(12);
    for (std::size_t i = 0; i < h.size(); i++)
        f << i << "  " << h[i] << "\n";
}

static Vec gauss(Mat A, Vec b)
{
    const int n = static_cast<int>(A.size());
    for (int k = 0; k < n; k++)
    {
        int p = k;
        for (int i = k + 1; i < n; i++)
            if (std::abs(A[i][k]) > std::abs(A[p][k])) p = i;

        std::swap(A[k], A[p]);
        std::swap(b[k], b[p]);

        if (std::abs(A[k][k]) < EPS)
            throw std::runtime_error("singular matrix");

        for (int i = k + 1; i < n; i++)
        {
            double m = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }

    Vec x(n);
    for (int i = n - 1; i >= 0; i--)
    {
        double s = b[i];
        for (int j = i + 1; j < n; j++) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return x;
}

static Vec lu_solve(const Mat& A, const Vec& b)
{
    const int n = static_cast<int>(A.size());
    Mat LU = A;
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);

    for (int k = 0; k < n; k++)
    {
        int p = k;
        for (int i = k + 1; i < n; i++)
            if (std::abs(LU[i][k]) > std::abs(LU[p][k])) p = i;

        std::swap(LU[k], LU[p]);
        std::swap(perm[k], perm[p]);

        if (std::abs(LU[k][k]) < EPS)
            throw std::runtime_error("singular matrix");

        for (int i = k + 1; i < n; i++)
        {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < n; j++)
                LU[i][j] -= LU[i][k] * LU[k][j];
        }
    }

    Vec y(n);
    for (int i = 0; i < n; i++)
    {
        double s = b[perm[i]];
        for (int j = 0; j < i; j++) s -= LU[i][j] * y[j];
        y[i] = s;
    }

    Vec x(n);
    for (int i = n - 1; i >= 0; i--)
    {
        double s = y[i];
        for (int j = i + 1; j < n; j++) s -= LU[i][j] * x[j];
        x[i] = s / LU[i][i];
    }
    return x;
}

static Vec jacobi(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0), xn(n);

    for (int it = 0; it <= maxIter; it++)
    {
        double r = residual_norm(A, x, b);
        hist.push_back(r);
        if (r < tol) break;

        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int j = 0; j < n; j++)
                if (j != i) s -= A[i][j] * x[j];
            xn[i] = s / A[i][i];
        }
        x.swap(xn);
    }
    return x;
}

static Vec seidel(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);

    for (int it = 0; it <= maxIter; it++)
    {
        double r = residual_norm(A, x, b);
        hist.push_back(r);
        if (r < tol) break;

        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int j = 0; j < n; j++)
                if (j != i) s -= A[i][j] * x[j];
            x[i] = s / A[i][i];
        }
    }
    return x;
}

static Vec sor(const Mat& A, const Vec& b, double omega, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);

    for (int it = 0; it <= maxIter; it++)
    {
        double r = residual_norm(A, x, b);
        hist.push_back(r);
        if (r < tol) break;

        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int j = 0; j < n; j++)
                if (j != i) s -= A[i][j] * x[j];

            double gs = s / A[i][i];
            x[i] = (1.0 - omega) * x[i] + omega * gs;
        }
    }
    return x;
}

static Vec gradient_descent(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);

    for (int it = 0; it <= maxIter; it++)
    {
        Vec r = sub(b, matvec(A, x));
        double rn = norm2(r);
        hist.push_back(rn);
        if (rn < tol) break;

        Vec Ar = matvec(A, r);
        double alpha = dot(r, r) / dot(Ar, r);
        for (int i = 0; i < n; i++) x[i] += alpha * r[i];
    }
    return x;
}

static Vec min_residual(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);

    for (int it = 0; it <= maxIter; it++)
    {
        Vec r = sub(b, matvec(A, x));
        double rn = norm2(r);
        hist.push_back(rn);
        if (rn < tol) break;

        Vec Ar = matvec(A, r);
        double alpha = dot(r, Ar) / dot(Ar, Ar);
        for (int i = 0; i < n; i++) x[i] += alpha * r[i];
    }
    return x;
}

static Vec conjugate_gradient(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);
    Vec r = sub(b, matvec(A, x));
    Vec p = r;
    double rr = dot(r, r);

    for (int it = 0; it <= maxIter; it++)
    {
        double rn = std::sqrt(rr);
        hist.push_back(rn);
        if (rn < tol) break;

        Vec Ap = matvec(A, p);
        double alpha = rr / dot(p, Ap);

        for (int i = 0; i < n; i++) x[i] += alpha * p[i];
        for (int i = 0; i < n; i++) r[i] -= alpha * Ap[i];

        double rr_new = dot(r, r);
        double beta = rr_new / rr;

        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rr = rr_new;
    }
    return x;
}

static Vec bicgstab(const Mat& A, const Vec& b, double tol, int maxIter, std::vector<double>& hist)
{
    const int n = static_cast<int>(A.size());
    Vec x(n, 0.0);
    Vec r = sub(b, matvec(A, x));
    Vec r0 = r;
    double rho_old = 1.0, alpha = 1.0, omega = 1.0;
    Vec v(n, 0.0), p(n, 0.0);

    double rn = norm2(r);
    hist.push_back(rn);
    if (rn < tol) return x;

    for (int it = 1; it <= maxIter; it++)
    {
        double rho = dot(r0, r);

        double beta = (rho / rho_old) * (alpha / omega);
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

        v = matvec(A, p);
        alpha = rho / dot(r0, v);

        Vec s(n);
        for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];

        double sn = norm2(s);
        if (sn < tol)
        {
            for (int i = 0; i < n; i++) x[i] += alpha * p[i];
            hist.push_back(sn);
            break;
        }

        Vec t = matvec(A, s);
        omega = dot(t, s) / dot(t, t);

        for (int i = 0; i < n; i++) x[i] += alpha * p[i] + omega * s[i];
        for (int i = 0; i < n; i++) r[i] = s[i] - omega * t[i];

        rho_old = rho;
        rn = norm2(r);
        hist.push_back(rn);
        if (rn < tol) break;
    }
    return x;
}

struct SLAE
{
    int n = 0;
    Mat A;
    Vec b;
};

static SLAE read_slae(const std::string& fname)
{
    std::ifstream f(fname);
    if (!f) throw std::runtime_error("cannot open file: " + fname);

    Vec nums;
    std::string line;
    while (std::getline(f, line))
    {
        auto pos = line.find('#');
        if (pos != std::string::npos) line.erase(pos);
        std::istringstream iss(line);
        double v;
        while (iss >> v) nums.push_back(v);
    }

    if (nums.empty())
        throw std::runtime_error("file is empty");

    std::size_t idx = 0;
    int n = static_cast<int>(nums[idx++]);
    if (n <= 0)
        throw std::runtime_error("invalid dimension n");

    const std::size_t need = static_cast<std::size_t>(n) * (n + 1);
    if (nums.size() < 1 + need)
        throw std::runtime_error("not enough numbers for A and b");

    SLAE s;
    s.n = n;
    s.A.assign(n, Vec(n, 0.0));
    s.b.assign(n, 0.0);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            s.A[i][j] = nums[idx++];
    for (int i = 0; i < n; ++i)
        s.b[i] = nums[idx++];
        
    return s;
}

static void print_vector(const std::string& name, const Vec& x)
{
    std::cout << name << " = [";
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (i) std::cout << ", ";
        std::cout << std::setprecision(8) << std::fixed << x[i];
    }
    std::cout << "]\n";
}

static void print_solution_block(const std::string& name, const Vec& x)
{
    std::cout << "\n--- " << name << " ---\n";
    print_vector("x", x);
}

int main(int argc, char** argv)
{
    std::string input_file;

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (a == "-i" || a == "--input")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "expected filename after " << a << "\n";
                return 1;
            }
            input_file = argv[++i];
        }
    }

    if (input_file.empty())
    {
        std::cerr << "no input file provided (use -i)\n";
        return 1;
    }

    SLAE s;
    try
    {
        s = read_slae(input_file);
        std::cout << "SLAE loaded from: " << input_file << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    const Mat& A = s.A;
    const Vec& b = s.b;
    const int n = s.n;
    const double tol = 1e-10;
    const int maxIter = 100000;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "n = " << n << ", tol = " << tol << "\n";

    auto report = [&](const std::string& name, const Vec& x)
    {
        double res = residual_norm(A, x, b);
        std::cout << std::scientific << std::setprecision(4) << "r = " << res << "\n";
        print_solution_block(name, x);
    };

    try { report("Gauss", gauss(A, b)); }
    catch (const std::exception& e) { std::cout << "Gauss: " << e.what() << "\n"; }

    try { report("LU", lu_solve(A, b)); }
    catch (const std::exception& e) { std::cout << "LU: " << e.what() << "\n"; }

    std::vector<double> hj, hs, hso, hg, hm, hcg, hb;

    {
        Vec x = jacobi(A, b, tol, maxIter, hj);
        report("Jacobi", x);
        std::cout << "   iterations: " << hj.size() - 1 << "\n";
    }
    {
        Vec x = seidel(A, b, tol, maxIter, hs);
        report("Seidel", x);
        std::cout << "   iterations: " << hs.size() - 1 << "\n";
    }
    {
        Vec x = sor(A, b, 1.5, tol, maxIter, hso);
        report("SOR (w=1.5)", x);
        std::cout << "   iterations: " << hso.size() - 1 << "\n";
    }
    {
        Vec x = gradient_descent(A, b, tol, maxIter, hg);
        report("Gradient descent", x);
        std::cout << "   iterations: " << hg.size() - 1 << "\n";
    }
    {
        Vec x = min_residual(A, b, tol, maxIter, hm);
        report("Minimal residual", x);
        std::cout << "   iterations: " << hm.size() - 1 << "\n";
    }
    {
        Vec x = conjugate_gradient(A, b, tol, maxIter, hcg);
        report("Conjugate gradient", x);
        std::cout << "   iterations: " << hcg.size() - 1 << "\n";
    }
    {
        Vec x = bicgstab(A, b, tol, maxIter, hb);
        report("BiCGStab", x);
        std::cout << "   iterations: " << hb.size() - 1 << "\n";
    }

    write_history("res_jacobi.txt", hj);
    write_history("res_seidel.txt", hs);
    write_history("res_sor.txt", hso);
    write_history("res_grad.txt", hg);
    write_history("res_minres.txt", hm);
    write_history("res_cg.txt", hcg);
    write_history("res_bicgstab.txt", hb);

    return 0;
}
