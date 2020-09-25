#include <iostream>
#include <vector>
#include <set>
#include <time.h>
#include <string>
#include <iomanip>
#include <chrono>
#include <bitset>
#include <array>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <fstream>


using namespace std;

const int k = 4;
const int m = 1383, n = 12, nn = 6;
int d = 0;
std::array<int, nn> MaxElem;
using TBitRow = std::bitset<n>;
using TBitCol = std::bitset<m>;

struct TMatrix
{
    vector<bitset<n * k>> row_mat = vector<bitset<n * k>>(m);
    vector<bitset<m>> col_mat = vector<bitset<m>>(n * k);
};

using TIntMatrix = std::array<std::array<int, n>, m>;


class ISolver {
protected:
    struct TContext {
        TContext(int size)
            : S(size)
            {}

        TBitCol D;
        TBitRow C, H;
        array<int, n> Sigma;
        vector<bitset<k - 1>> S;
    } Context;

    vector<int> Count;
    TIntMatrix Mat;

public:
    ISolver(TIntMatrix mat)
        : Context(n)
        , Mat(mat)
    {
        for (size_t i = 0; i < m; ++i) {
            if (i < d) {
                Context.D.set(i);
            } else {
                Context.D.reset(i);
            }
            Count.push_back(accumulate(Mat[i].begin(), Mat[i].end(), 0));
        }
        for (size_t j = 0; j < n; ++j) {
            for (size_t p = 0; p < k - 1; ++p) {
                Context.S[j].set(p);
            }
            
            Context.C.set(j);
        }
    }

    virtual std::pair<vector<TBitRow>, vector<array<int, n>>> solve() = 0;
};

class SmartSolver : public ISolver {
    size_t GetWithMinWeight() {
        ssize_t minIdx = -1;
        for (size_t idx = 0; idx < m; ++idx) {
            if (Context.D[idx]) {
                if (minIdx == -1) {
                    minIdx = idx;
                } else if (Count[idx] < Count[minIdx]) {
                    minIdx = idx;
                }
            }
        }
        return minIdx;
    }
    
    void BuildSubtreeRUNC(vector<TBitRow>& ans, vector<array<int, n>>& ans_sigma, int depth = 0) {

        size_t i = GetWithMinWeight();
        //cout << "min " << i << ' ' << depth << endl;
        for (int colIdx = 0; colIdx < n; ++colIdx) {

            if (!Context.C[colIdx] || Mat[i][colIdx] == 0) {
                continue;
            }
            bool good = false;
            for (int t = Mat[i][colIdx] - 1; t >= 0; --t) {
                if (t > k - 2) {
                    continue;
                }
                if (Context.S[colIdx].test(t)) {
                    good = true;
                    break;
                }
                
            }
            if (!good) {
                continue;
            }
            Context.C.reset(colIdx); // (6)
            TContext savedContext(Context);
            
            for (int x = Mat[i][colIdx] - 1; x >= 0; --x) {
                if (colIdx - nn >= 0 && Context.H.test(colIdx - nn) && Context.Sigma[colIdx - nn] < MaxElem[colIdx - nn] - x) {
                    continue;
                }
                if (colIdx + nn < n && Context.H.test(colIdx + nn) && MaxElem[colIdx] - Context.Sigma[colIdx + nn] > x) {

                    continue;
                }
                if (x <= k - 2 && Context.S[colIdx].test(x)) {
                    Context.S[colIdx].reset(x);

                    Context.H.set(colIdx);
                    Context.Sigma[colIdx] = x;

                    for (int row = 0; row < m; ++row) {
                        if (Context.D[row] && Mat[row][colIdx] > x) {
                            Context.D.reset(row);
                        }
                    }
                    if (Context.D.count() == 0) {
                        ans.push_back(std::move(Context.H));
                        ans_sigma.push_back(std::move(Context.Sigma));
                        Context = std::move(savedContext);
                        Context.S[colIdx].reset(x);
                        break;
                    } else {
                        for (int p = 0; p < n; ++p) {
                            if (!Context.C[p]) {
                                continue;
                            }
                            for (int v = 0; v <= k - 2; ++v) {
                                if (Context.S[p].test(v)) {
                                    Context.Sigma[p] = v;
                                    Context.H.set(p); // (*)
                                    bool ok = CheckSigma();
                                    if (!ok) {
                                        Context.S[p].reset(v);
                                    }
                                    Context.H.reset(p); // rollback changes from (*)
                                }
                            }
                        }
                        BuildSubtreeRUNC(ans, ans_sigma, depth + 1);
                    }
                    Context = savedContext;
                    Context.S[colIdx].reset(x);// no move allowed
                }
            }
            Context.C.set(colIdx); // rollback (7)
        }
    }

    bool CheckSigma() {
        for (int col = 0; col < n; ++col) {
            if (!Context.H[col]) {
                continue;
            }
            bool found_correct_row = false;
            for (int row = 0; row < m; ++row) {
                bool correct = true;
                for (int idx = 0; idx < n; ++idx) {
                    if (!Context.H[idx]) {
                        continue;
                    }
                    if (idx == col) {
                        if (Mat[row][idx] != (1 + Context.Sigma[idx])) {
                            correct = false;
                            break;
                        }
                    } else if (Mat[row][idx] > Context.Sigma[idx]) {
                        correct = false;
                        break;
                    }
                }
                if (correct) {
                    found_correct_row = true;
                    break;
                }
            }
            if (!found_correct_row) {
                return false;
            }
        }
        return true;
    }

public:
    SmartSolver(TIntMatrix mat) : ISolver(mat) {}

    std::pair<vector<TBitRow>, vector<array<int, n>>> solve() override {
        vector<TBitRow> ans;
        vector<array<int, n>> ans_sigma;
        BuildSubtreeRUNC(ans, ans_sigma);
        return std::make_pair(ans, ans_sigma);
    }
};



vector<TIntMatrix> GenerateIntMatrixes() {
    int count = 10;
    vector<TIntMatrix> res;
    while (count--) {
        res.emplace_back();
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res.back()[i][j] = rand() % k;
            }
        }
    }
    return res;
}

TIntMatrix ReadIntMatrixes() {

    TIntMatrix res;
    std::string line;
    std::ifstream in("in.txt");
    int i = 0;
    if (in.is_open())
    {
        int num;
        for (int j = 0; j < nn; ++j) {
            in >> num;
            MaxElem[j] = num;
        }
        while (!in.eof())
        {
            in >> num;
            res[i / n][i % n] = num;
            ++i;
        }
        
    }
    in.close();
    d = i / n;
    return res;
}


TMatrix Convert(const TIntMatrix& input) {
    TMatrix res;
    for (size_t i = 0; i < m; ++i) {
        size_t bitIdx = 0;
        for (size_t j = 0; j < n; ++j) {
            for (size_t r = 0; r < k; ++r) {
                res.row_mat[i].set(bitIdx) = (r <= input.at(i).at(j));
                res.col_mat[bitIdx].set(i) = (r <= input.at(i).at(j));
                ++bitIdx;
            }
        }
    }
    return res;
}

vector<TMatrix> Convert(const vector<TIntMatrix>& input) {
    vector<TMatrix> res(input.size());
    for (size_t idx = 0; idx < res.size(); ++idx) {
        res[idx] = Convert(input[idx]);
    }
    return res;
}

string GenerateStringFromCoverages(const vector<TBitRow>& ans) {
    set<string> all;
    for (const auto& coverage : ans) {
        string coverageStr;
        for (int x = 0; x < n; ++x) {
            if (coverage[x]) {
                coverageStr += (to_string(x) + " ");
            }
        }
        all.insert(std::move(coverageStr));
    }
    string result;
    for (auto coverageStr : all) {
        result += "#" + coverageStr;
    }
    return result;
}

int main() {
    const auto& mat = ReadIntMatrixes();
    
   /* if (n * k > 20){
        std::cout << 'No correct test';
    }
    for (const auto& mat : input) {
        SmartSolver smart(mat);
        StupidSolver stupid(mat);
        auto kek1 = smart.solve();
        auto kek2 = stupid.solve();
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << mat[i][j] << ' ';
            }
            cout << endl;
        }
        if (GenerateStringFromCoverages(kek1) != GenerateStringFromCoverages(kek2)) {
            cout << "BAD" << endl;
        }
        //cout << "smart - " << GenerateStringFromCoverages(kek1) << endl;
        //cout << "stupid - " << GenerateStringFromCoverages(kek2) << endl;
    }*/
    std::chrono::nanoseconds total;
    SmartSolver smart(mat);
    auto start = std::chrono::system_clock::now();
    auto answer_pair = smart.solve();
    auto end = std::chrono::system_clock::now();
    total += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    auto answer = answer_pair.first;
    auto answer_sigma = answer_pair.second;
    cout << "thre";
    std::ofstream out;
    out.open("out.csv");
    if (out.is_open())
    {
        for  (const auto& coverage : answer) {
            for (int col = 0; col < n; ++col) {
                if (coverage[col])
                    out << 1;
                if (!coverage[col])
                    out << 0;
                if (col < n - 1)
                    out << ',';
                    
                
            }
            out << std::endl;
        }
    }
    out.close();
    std::ofstream out1;
    out1.open("out_sigma.csv");
    if (out1.is_open())
    {
        for  (const auto& sigma : answer_sigma) {
            for (int col = 0; col < n; ++col) {
                out1 << sigma[col];
                if (col < n - 1)
                    out1 << ',';
            }
            out1 << std::endl;
        }
    }
    out1.close();
    std::cout << "End of program" << std::endl;


    auto milliseconds = total / 1'000'000;

    cout << "m " << m << " n " << n << endl;
    cout << setprecision(10)
        << total.count() / 1'000'000'000. << endl
        << milliseconds.count() / 1'000.f << endl;
	
}
