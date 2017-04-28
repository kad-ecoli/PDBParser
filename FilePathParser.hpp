/* parse path of file */
#include<string>
#include<cstring>

/* extract the basename from a file path */
string basename_no_ext(const char *pdb_file,bool suppress_ext=true)
{
    string filename;
    int startindex=0;
    int endindex=strlen(pdb_file);
    for (int i=0;i<strlen(pdb_file);i++){
        if (pdb_file[i]=='/'||pdb_file[i]=='\\'){
            startindex=i+1;
            endindex=strlen(pdb_file);
        }
        if (pdb_file[i]=='.' && endindex==strlen(pdb_file) && suppress_ext){
            endindex=i;
        }
    }   
    for (int i=startindex;i<endindex;i++) filename+=pdb_file[i];
    return filename;
}

/* extract the full filename minus file extension */
string filename_no_ext(const char *pdb_file)
{
    string filename;
    int endindex=strlen(pdb_file);
    for (int i=endindex-1;i>=0;i--){
        if (pdb_file[i]=='/'||pdb_file[i]=='\\') break;
        if (pdb_file[i]=='.') endindex=i;
    }   
    for (int i=0;i<endindex;i++) filename+=pdb_file[i];
    return filename;
}
