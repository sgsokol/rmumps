/unused variable/ {
   #print "unused var";
   split($1, a, ":");
   ("find -name " a[1]) | getline path;
   print path, a[2], $5;
   cmd_sed="sed -i.sedbck '" a[2] "{s/\\s*$//; s:\\(.*[^/]\\)\\(/[*].*[*]/\\)\\?:/*\\1*/\\2:}' " path;
   print cmd_sed;
   (cmd_sed) | getline status;
   ("geany " path " +" a[2]) | getline status;
   #exit 0;
}
/__escape__ set but not used \[-Wunused-but-set-variable\]/ {
   #print "set but unused"
   split($1, a, ":");
   ("find -name " a[1]) | getline path;
   print path, a[2], $4;
   #("geany " path " +" a[2]) | getline status;
}
