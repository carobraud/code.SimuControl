
c ******** POUR ERREUR DANS LECTURE NETCDF
      subroutine handle_err(status)
c      include 'netcdf.inc'

      integer status

      if (status .ne. nf_noerr) then
        print *, nf_strerror(status)
        stop 'Problem with NetCDF'
      endif

      return
      end

