c ******** Subroutine ECRITURE NETCDF (uniquement plan PIV 3C)

      subroutine Write_Ncdf(pathresult,fileout,nx,ny,nz,x,y,z,vx,vy,vz)
      integer idcount(3)
      integer start(3),count(3)
      dimension x(*),y(*),z(*)
      dimension vx(*),vy(*),vz(*)      
      character pathresult*100,fileout*100
c     ncase=1 : champs moyen
c
c     Utiliser fonctions définies dans netcdf.inc
c     1- Create netCDF file : NF_CREATE (!!!ATTENTION!!! NF_CLOBBER écrase le fichier du meme nom, mettre NF_NOCLOBBER sinon) 
c     2- Define dimension for output : NF_DEF_DIM - one per dimension
c     3- Define dimension Variable : NF_DEF_VAR - one per dimension
c     4- Define data Variable : NF_DEF_VAR - one per dimension
c     5- End the 'define' mode in netCDF : NF_ENDDEF
c     6- Define the 'shape' of the data to be output.
c     7- Output the values of the dimension variables : NF_PUT_VARA_DOUBLE
c     8- Output the values of the data variables : NF_PUT_VARA_DOUBLE
c     9- Close the netCDF file : NF_CLOSE

c     !!! ATTENTION les IDs doivent etre lus dans l'ordre (inverse de celui donne par ncdump: x,y,z,t)

c     1- OPEN netCDF file

      include "netcdf.inc"
      istatus=nf_create(pathresult(1:len_trim(pathresult))//fileout(1:
     &	len_trim(fileout)),nf_CLOBBER,ncid)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'CREATE FILE :',pathresult(1:len_trim(pathresult))//fileout
     &     (1:len_trim(fileout))

c     2- Define dimension for output
      istatus=nf_def_dim(ncid,'dim_x',nx,idimx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'Id dim x',nx
      istatus=nf_def_dim(ncid,'dim_y',ny,idimy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'Id dim y',ny
      istatus=nf_def_dim(ncid,'dim_z',nz,idimz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'Id dim z',nz

c     3- Define Dimension Variable
      istatus=nf_def_var(ncid,'grid_x',nf_double,1,idimx,ivaridx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'grid x',ivaridx
      istatus=nf_def_var(ncid,'grid_y',nf_double,1,idimy,ivaridy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'grid y',ivaridy
      istatus=nf_def_var(ncid,'grid_z',nf_double,1,idimz,ivaridz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'grid z',ivaridz

c     4- Define Data Variable
      idstart=3
      idcount(1)=idimx
      idcount(2)=idimy
      idcount(3)=idimz

      istatus=nf_def_var(ncid,'vel_x',nf_double,idstart,idcount,ivar
     &	idvx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vel x',ivaridvx
      istatus=nf_def_var(ncid,'vel_y',nf_double,idstart,idcount,ivar
     & 	idvy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vel y',ivaridvy
      istatus=nf_def_var(ncid,'vel_z',nf_double,idstart,idcount,ivar
     &	idvz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vel z',ivaridvz

c     5- end the 'define' mode in netCDF
      istatus=nf_ENDDEF (ncid)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

c     6/7- Output the values of the dimension variables : 

c      istatus=nf_put_vara_double(ncid,ivaridx,1,nx,x)
c      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'x',istatus
c      istatus=nf_put_vara_double(ncid,ivaridy,1,ny,y)
c      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'y',istatus
c      istatus=nf_put_vara_double(ncid,ivaridz,1,nz,z)
c      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'z',istatus

c     8- Output the values of the data variables : 
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=nx
      count(2)=ny
      count(3)=nz
      istatus=nf_put_vara_double(ncid,ivaridvx,start,count,vx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      istatus=nf_put_vara_double(ncid,ivaridvy,start,count,vy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      istatus=nf_put_vara_double(ncid,ivaridvz,start,count,vz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

c     9- Close the file up
 111  istatus=nf_close(ncid)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

      return
      end

      
