c
c     CODE WHICH PRODUCES A FRACTION OF CONTROL.DAT
c    
      A=0.020
      N=50
      DN= 0.02
      open (1, file='filenames.dat')
      write (1,99) N 
      DO i=1,N
        write (1,100) A  
         A=A + DN 
      ENDDO
 99   format (I6)
 100  format (F6.3)
      END
