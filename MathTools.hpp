/* functions for numpy style matrix manipulation */
#if !defined(MATHTOOLS_HPP)
#define MATHTOOLS_HPP 1

#include <iostream> 
#include <iomanip>
#include <vector>

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

using namespace std;

/* print the content of vector 
 * max_num is the max number of element to print */
void print_vector(const vector<int>& vec,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int vec_len=vec.size();
    cout<<"[";
    for(int i=0;i<vec_len;i++)
    {
        if (i==(vec_len-1)) cout<<vec[i]<<"]"<<endl;
        else if (i<max_num/2 || vec_len-i <=max_num/2) cout<<vec[i]<<", ";
        else if (i==max_num/2 && vec_len>max_num) cout<<" ... ";
    }
    if (vec_len==0) cout<<"]"<<endl;
}

/* print the content of vector */
void print_vector(const vector<double>& vec,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int vec_len=vec.size();
    cout<<"[";
    for(int i=0;i<vec_len;i++)
    {
        if (i==(vec_len-1)) cout<<setprecision(4)<<vec[i]<<"]"<<endl;
        else if (i<max_num/2 || vec_len-i <=max_num/2) cout<<setprecision(4)<<vec[i]<<", ";
        else if (i==max_num/2 && vec_len>max_num) cout<<" ... ";
    }
    if (vec_len==0) cout<<"]"<<endl;
}

/* print the content of vector */
void print_vector(const vector<string>& vec,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int vec_len=vec.size();
    cout<<"[";
    for(int i=0;i<vec_len;i++)
    {
        if (i==(vec_len-1)) cout<<vec[i]<<"]"<<endl;
        else if (i<max_num/2 || vec_len-i <=max_num/2) cout<<vec[i]<<", ";
        else if (i==max_num/2 && vec_len>max_num) cout<<" ... ";
    }
    if (vec_len==0) cout<<"]"<<endl;
}

/* print the content of matrix 
 * max_num is the maximum number of row/columns to print */
void print_matrix(const vector<vector<double> >& mat,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int m=mat.size();
    cout<<"["<<endl;
    for(int i=0;i<m;i++)
    {
        if (i<max_num/2 || m-i <=max_num/2) print_vector(mat[i],max_num);
        else if (i==max_num/2 && m>max_num) cout<<" ... "<<endl;
    }
    cout<<"]"<<endl;
}

/* print the content of matrix */
void print_matrix(const vector<vector<int> >& mat,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int m=mat.size();
    cout<<"["<<endl;
    for(int i=0;i<m;i++)
    {
        if (i<max_num/2 || m-i <=max_num/2) print_vector(mat[i],max_num);
        else if (i==max_num/2 && m>max_num) cout<<" ... "<<endl;
    }
    cout<<"]"<<endl;
}

/* print the content of matrix */
void print_matrix(const vector<vector<string> >& mat,int max_num=20)
{
    if (max_num % 2 !=0) max_num--; // must be even number
    int m=mat.size();
    cout<<"["<<endl;
    for(int i=0;i<m;i++)
    {
        if (i<max_num/2 || m-i <=max_num/2) print_vector(mat[i],max_num);
        else if (i==max_num/2 && m>max_num) cout<<" ... "<<endl;
    }
    cout<<"]"<<endl;
}

/* append or generate a list of interger number from start_num till end_num-1 */
void range(vector<int>& vec, const int start_num, const int end_num)
{
    for(int i=start_num;i<=end_num;i++) vec.push_back(i);
}
#endif
