#include <ctype.h>
#include <string.h>
#include <sys/utsname.h>
#define NUL  0
#define DOT 46
#define MAX 63

int ihost_(host_name) /* DEC, SGI, Sun, Linux */
char       host_name[MAX] ;
{
  int    i = 0          ;
  struct utsname opname ;
  if (uname (&opname) == -1)
  {
    return (-1)         ;
  }
  else
  {
    if (strcpy (host_name, opname.nodename) == NULL)
    {
      return (-2) ;
    }
    else
    {
      host_name[i] = toupper (host_name[i])   ;
      while ( host_name[i] != NUL && host_name[i] != DOT && i < MAX )
      {
        host_name[i] = toupper (host_name[i]) ;
        i++                                   ;
      }
      host_name[i] = NUL ;
      return (i)         ;
    }
  }
}
