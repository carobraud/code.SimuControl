
c ******** Subroutine ECRITURE NETCDF (uniquement plan PIV 3C)

      subroutine Write_Ncdf(pathresult,fileout,nx,ny,x,y,vx,vy,wz,mask)
      integer idcount(2)
      integer start(2),count(2)
      dimension x(*),y(*)
      dimension vx(*),vy(*),wz(*),mask(*)      
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

c     3- Define Dimension Variable
      istatus=nf_def_var(ncid,'grid_x',nf_real,1,idimx,ivaridx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'grid x',ivaridx
      istatus=nf_def_var(ncid,'grid_y',nf_real,1,idimy,ivaridy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'grid y',ivaridy

c     4- Define Data Variable
      idstart=2
      idcount(1)=idimx
      idcount(2)=idimy

      istatus=nf_def_var(ncid,'vel_x',nf_real,idstart,idcount,ivar
     &	idvx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vel x',ivaridvx
      istatus=nf_def_var(ncid,'vel_y',nf_real,idstart,idcount,ivar
     & 	idvy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vel y',ivaridvy
      istatus=nf_def_var(ncid,'vort_z',nf_real,idstart,idcount,ivar
     & 	idwz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'vort z',ivaridwz
        istatus=nf_def_var(ncid,'mask',nf_real,idstart,idcount,ivar
     &     idmask)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      print*,'mask',ivaridmask
      
c     5- end the 'define' mode in netCDF
      istatus=nf_ENDDEF (ncid)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

c     6/7- Output the values of the dimension variables : 

      istatus=nf_put_vara_real(ncid,ivaridx,1,nx,x)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'x',istatus
      istatus=nf_put_vara_real(ncid,ivaridy,1,ny,y)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'y',istatus
c      istatus=nf_put_vara_double(ncid,ivaridz,1,nz,z)
c      if (istatus.ne.nf_noerr) call handle_err(istatus)
c      print*,'z',istatus

c     8- Output the values of the data variables : 
      start(1)=1
      start(2)=1
      count(1)=nx
      count(2)=ny
      istatus=nf_put_vara_real(ncid,ivaridvx,start,count,vx)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      istatus=nf_put_vara_real(ncid,ivaridvy,start,count,vy)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      istatus=nf_put_vara_real(ncid,ivaridwz,start,count,wz)
      if (istatus.ne.nf_noerr) call handle_err(istatus)
      istatus=nf_put_vara_real(ncid,ivaridmask,start,count,mask)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

c     9- Close the file up
 111  istatus=nf_close(ncid)
      if (istatus.ne.nf_noerr) call handle_err(istatus)

      return
      end

      
