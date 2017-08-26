/* functions for python style string manipulation */
#if !defined(STRINGTOOLS_HPP)
#define STRINGTOOLS_HPP 1

#include <string> 
#include <vector>

using namespace std;

/* join string_vec into a long string using seperator "sep" */
string join(const string sep, const vector<string>& string_vec)
{
    if (string_vec.size()==0) return "";
    string joined_str=string_vec[0];
    for (int s=1;s<string_vec.size();s++) joined_str+=sep+string_vec[s];
    return joined_str;
}
#endif
