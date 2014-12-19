#ifndef DO_MKDIR_H
#define DO_MKDIR_H

#include <errno.h> // global 'errno' parameter
#include <sys/stat.h> // mkdir, stat

// do_mkdir() by Jonathan Leffler, Public Domain
// http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux
int do_mkdir(const char *path, mode_t mode)
{
  struct stat st;
  int status = 0;

  if (stat(path, &st) != 0)
    {
      // Directory does not exist. EEXIST for race condition
      if (mkdir(path, mode) != 0 && errno != EEXIST)
        status = -1;
    }
  else if (!S_ISDIR(st.st_mode))
    {
      errno = ENOTDIR;
      status = -1;
    }

  return status;
}

#endif
