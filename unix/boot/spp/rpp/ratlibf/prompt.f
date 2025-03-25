      subroutine prompt (str, buf, fd)
      integer str(100), buf(100)
      integer fd
      logical isatty
      if (.not.(isatty(fd)))goto 23000
      call putlin (str, fd)
      call rfflus(fd)
23000 continue
      call getlin (buf, fd)
      return
      end
