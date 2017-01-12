        Program  TJ_line
		    

          implicit none

          integer::i,j,k,n,npts1,npts2,npts3
          real::TJ,x,y,z
          integer,dimension(128,128,1)::TJL


          open(5,file="2d_GB_list_test.txt",status="old",action="read")
          open(1,file="output.vtk",status="replace",action="write")

          Do n=1,1292

             read(5,*) z,y,x,TJ

             i=int(x)+1
             j=int(y)+1
!             k=int(z)
!             write(*,*) i,j
             TJL(i,j,1)=int(TJ)

          end Do



          npts1=128
          npts2=128
          npts3=1


         write(1,'(a)') '# vtk DataFile Version 2.0'
         write(1,'(a)') 'Equivalent strain field'
         write(1,'(a)') 'ASCII'
         write(1,'(a)') 'DATASET STRUCTURED_POINTS'
         write(1,'(a,3(1x,i5))') 'DIMENSIONS',npts1+1,npts2+1,npts3+1
         write(1,'(a,3(1x,f10.5))') 'ASPECT_RATIO',1.0,1.0,1.0
         write(1,'(a)') 'ORIGIN 0 0 0'
         write(1,'(a,1x,i10)') 'CELL_DATA',npts1*npts2*npts3
         write(1,'(a)') 'SCALARS phase Int'
         write(1,'(a)') 'LOOKUP_TABLE default'

         Do k=1,npts3
          Do j=1,npts2
            Do i=1,npts1

               write(1,*) TJL(i,j,k)
              
             enddo
          enddo
       end Do

     end Program TJ_line
     

          
          
          
