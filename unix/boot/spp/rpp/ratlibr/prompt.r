include	defs

# prompt - write to/read from teletype

   subroutine prompt (str, buf, fd)
   character str(ARB), buf(ARB)
   filedes fd

   logical isatty

   if (isatty(fd))
	 {
	 call putlin (str, fd)
	 call flush (fd)
	 }
   call getlin (buf, fd)

   return
   end
