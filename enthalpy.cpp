// Enthalpy calculator, written by Albert Nikolay Bulik for CISC119, fall semester 2021. 

#include <vector>
#include <iostream>

namespace anb
{
    class Matrix
    {
    public:
        std::vector<std::vector<double>> row;

    public:
        Matrix(int rowSize, int colSize)
        {
            row = std::vector<std::vector<double>>(rowSize, std::vector<double>(colSize, 0));
        }

        Matrix(std::vector<std::vector<double>> data)
        {
            row = data;
        }
        
        void ReduceToRowEchelonForm()
        {
            int lead = 0;
            int rowCount = row.size();
            int columnCount = row[0].size();

            for (int r = 0; r < rowCount; r++)
            {
                if (columnCount <= lead) break;
                int i = r;

                while (row[i][lead] == 0) {
                    i++;

                    if (i == rowCount) 
                    {
                        i = r;
                        lead++;

                        if (columnCount == lead) {
                            lead--;
                            break;
                        }
                    }
                }

                for (int j = 0; j < columnCount; j++) 
                {
                    int temp = row[r][j];
                    row[r][j] = row[i][j];
                    row[i][j] = temp;
                }

                int div = row[r][lead];
                if(div != 0) 
                {
                    for (int j = 0; j < columnCount; j++)
                    {
                        row[r][j] /= div;
                    }
                }

                for (int j = 0; j < rowCount; j++)
                {
                    if (j != r)
                    {
                        int sub = row[j][lead];
                        for (int k = 0; k < columnCount; k++) row[j][k] -= (sub * row[r][k]);
                    }
                }
                lead++;
            }
        } // end row reduction function
    }; // end class 
}

// overload the << operator when the object to the right is type matrix
std::ostream &operator<<(std::ostream &os, anb::Matrix const &m)
{
    os << "(Matrix dimensions: " << m.row.size() << " x " << m.row[0].size() << ")" << std::endl;
    for (auto col: m.row)
    {
        for (auto entry: col)
        {
            os << entry << "\t";
        }
        os << "\n";
    }
    os << "(End matrix output)" << std::endl;
    return os;
};

// unused function to test matrix row reduction
static void TestMatrix() 
{

    {
        auto mat = anb::Matrix({
            {  1, 2, -1,  -4 },
            { 2, 3, -1, -11 },
            { -2, 0, -3,  22 }
        });
        mat.ReduceToRowEchelonForm();
        std::cout << mat << std::endl;
    }

    {
        auto mat = anb::Matrix(3, 4);
        mat.row[0] = {  1, 2, -1,  -4 };
        mat.row[1] = { 2, 3, -1, -11 };
        mat.row[2] = { -2, 0, -3,  22 };
        mat.ReduceToRowEchelonForm();
        std::cout << mat << std::endl;
    }

}
// end Matrix class

// #include <iostream>
// #include <vector>
#include <algorithm>
#include <string>

using namespace std;

// function to split string
vector<string> split(const string& s)
{
    vector<string> ret;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) 
    {
        // ignore leading blanks
        while (i != s.size() && isspace(s[i]))
            ++i; 
        
        // find end of next word
        string_size j = i; 
        while (j != s.size() && !isspace(s[j]))
            ++j;
        
        // if we found some non-whitespace characters
        if (i != j) 
        {
            // copy from s starting at i and taking j - i chars
            ret.push_back(s.substr(i, j - i));
            i = j;
        }
    }

    return ret;
}

int main()
{
    cout    << "Please enter target reaction or 'help':\n\n" 
            << "E.g. P4O10 + 6H2O -> 4H3PO4\n" << endl; 

    // read and split first line of input
    string s;
    string::size_type h;
    getline(cin, s); 
    // this block is meant to remove lines of input that are purely whitespace
    {
        h = 0;
        // ignore leading blanks
        while (h != s.size() && isspace(s[h]))
            ++h;
        // cut off leading whitespace
        s = s.substr(h, s.size() - h);
        if (s.size() == 0) 
        {
            cout << "Invalid input." << endl;
            return 0;
        }
    }

    vector<string> targetRxn = split(s);

    // parse for 'help'
    s = targetRxn[0];
    transform(s.begin(), s.end(), s.begin(), ::toupper);
    if (!s.compare("HELP")) 
    {
        // run this if the first word is help
        cout << "Help text placeholder." << endl; 
    } 
    else // proceed with execution if 'help' is not mentioned
    {
        vector<string> formula; 
        vector<int> targetRxnCoefficients;
        bool arrowFlag = false; 

        // process all substrings of targetRxn
        for (vector<string>::size_type i = 0; i != targetRxn.size(); ++i) 
        {
            if (!targetRxn[i].compare("+")) continue;

            // check for arrow (->) in reaction equation 
            if (targetRxn[i].back() == '>') 
            {
                arrowFlag = true; 
            }
            else 
            {
                // determine coefficient of chemical formula
                string t;
                int coefficient; 
                if (isdigit(targetRxn[i].front())) 
                {
                    coefficient = stoi(targetRxn[i]);
                    // remove the coefficient in front of the chemical formula 
                    t = to_string(coefficient);
                    t = targetRxn[i].substr(t.size(), targetRxn[i].size() - t.size());
                }
                else 
                {
                    coefficient = 1; 
                    t = targetRxn[i]; 
                }
                formula.push_back(t);

                if (arrowFlag) coefficient *= -1;
                targetRxnCoefficients.push_back(coefficient);
            }
        } // end for

        cout    << "Please enter additional reactions, followed by enthalpy in kilojoules,\n" 
		        << "\tUse space and semicolon to separate reaction and enthalpy,\n"
		        << "\tCtrl-Z in Windows/Ctrl-D in Unix to end input:\n\n"
                << "E.g.\tP4 + 5O2 -> P4O10 ; -2984\n"
		        << "\t2H2 + O2 -> 2H2O ; -571.6\n"
		        << "\tP4 + 6H2 + 8O2 -> 4H3PO4 ; -2568.8\n\t(ctrl-z/ctrl-d)\n" << endl;
        
        vector<vector<string>> additionalRxns;
        // user inputs additional reactions
        while (getline(cin, s))
        {
            // this block is meant to remove lines of input that are purely whitespace
            {
                h = 0;
                // ignore leading blanks
                while (h != s.size() && isspace(s[h]))
                    ++h;
                // cut off leading whitespace
                s = s.substr(h, s.size() - h);
                if (s.size() == 0) continue;
            }

            vector<string> rxn = split(s);
            additionalRxns.push_back(rxn);
            bool enthalpyFlag = false; 
            // for each string-type entry in rxn[]
            for (vector<string>::size_type i = 0; i != rxn.size(); ++i)
            {
                if (!rxn[i].compare("+")) continue;
                if (rxn[i].back() == '>') continue;
                // ignore all values to the right of ';'
                if (!rxn[i].compare(";")) enthalpyFlag = true; 
                if (enthalpyFlag == true) continue;

                string t; 
                int coefficient; 
                if (isdigit(rxn[i].front())) 
                {
                    coefficient = stoi(rxn[i]);
                    // remove the coefficient in front of the chemical formula 
                    t = to_string(coefficient);
                    t = rxn[i].substr(t.size(), rxn[i].size() - t.size());
                }
                else t = rxn[i];
                // invariant: string t has the chemical formula only (no coefficient)

                bool isInSet = false; 
                // for each string-type entry in formula[]
                for (vector<string>::size_type j = 0; j != formula.size(); ++j)
                {
                    if (!t.compare(formula[j])) isInSet = true; 
                }
                if (!isInSet) formula.push_back(t); 
            }
        } // end while, bool enthalpyFlag is destroyed and reinitialized to false
        // invariant:   formula[] contains all unique chemical formulas
        //              additionalRxns[][] contains string vectors of reactions (except target rxn)

        // construct matrix
        auto coefficientMatrix = anb::Matrix(targetRxnCoefficients.size(), 0);
        for (int i = 0; i != targetRxnCoefficients.size(); ++i)
        {
            for (int j = 0; j != formula.size(); ++j)
            {
                bool hasVariable = false; 
                bool arrowFlag = false;
                for (int k = 0; k != additionalRxns[i].size(); ++k)
                {
                    if (!additionalRxns[i][k].compare("+")) continue;
                    if (!additionalRxns[i][k].compare(";")) break;
                    if (additionalRxns[i][k].back() == '>') 
                    {
                        arrowFlag = true; 
                        continue;
                    }

                    string t; 
                    int coefficient; 
                    if (isdigit(additionalRxns[i][k].front()))
                    {
                        coefficient = stoi(additionalRxns[i][k]);
                        // remove the coefficient in front of the chemical formula 
                        t = to_string(coefficient);
                        t = additionalRxns[i][k].substr(t.size(), additionalRxns[i][k].size() - t.size());
                    }
                    else 
                    {
                        coefficient = 1;
                        t = additionalRxns[i][k];
                    }
                    if (arrowFlag) coefficient *= -1;

                    if (!formula[j].compare(t)) 
                    {
                        coefficientMatrix.row[i].push_back(coefficient);
                        hasVariable = true;
                        break;
                    }
                } // inner
                if (!hasVariable)
                {
                    coefficientMatrix.row[i].push_back(0);
                }
            } // middle
            coefficientMatrix.row[i].push_back(stod(additionalRxns[i].back()));
        } // outer

        cout << "\n"; 
        // print column markers
        for (int i = 0; i != formula.size(); ++i) 
            cout << formula[i] << "\t";
        cout << "dH\n" << endl; 

        // print matrix before & after row reduction
        cout << coefficientMatrix << "\nRow reducing...\n";
        coefficientMatrix.ReduceToRowEchelonForm();
        cout << coefficientMatrix << endl;

        // multiply each row by corresponding values in targetRxnCoefficients
        for (int i = 0; i != targetRxnCoefficients.size(); ++i)
        {
            for (int j = 0; j != coefficientMatrix.row[i].size(); ++j) 
            {
                coefficientMatrix.row[i][j] *= targetRxnCoefficients[i]; 
            }
        } 
        cout << "Multiplying...\n" << coefficientMatrix << endl;

        // print column markers
        for (int i = 0; i != targetRxnCoefficients.size(); ++i)
            cout << "\t";
        for (int i = targetRxnCoefficients.size(); i != formula.size(); ++i) 
            cout << formula[i] << "\t";
        cout << "dH" << endl;

        for (int i = 0; i != targetRxnCoefficients.size(); ++i)
            cout << "\t";
        // sum vars on right
        bool balancedRxn = true; 
        for (int i = targetRxnCoefficients.size(); i != formula.size(); ++i)
        {
            int sum = 0;
            for (int j = 0; j != targetRxnCoefficients.size(); ++j)
                sum += coefficientMatrix.row[j][i];
            cout << sum << "\t";

            if (sum != 0) balancedRxn = false; 
        }

        // sum enthalpy
        double dH = 0.0; 
        for (int i = 0; i != targetRxnCoefficients.size(); ++i)
            dH += coefficientMatrix.row[i].back();
        cout    << dH << "\n\nThe enthalpy of this reaction is: " 
                << dH << " kJ\n";
        
        if (!balancedRxn)
        {
            cout    << "There are leftover reactants/products; " 
                    << "please provide additional equations." << endl;
        }
    } // end else

    // pause
    cout << "Press enter to exit." << endl;
    getchar();

    return 0; 
}