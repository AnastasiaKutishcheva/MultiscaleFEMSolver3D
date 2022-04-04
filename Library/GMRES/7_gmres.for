*****************************************************
	SUBROUTINE GMRS(IPTR,JPTR,LTR,UTR,DIAG,F,X,N,M,EL_NUMBER,EPS)
	IMPLICIT NONE
	!����������� ������, ����� �������� ��������
	!��������������� ������� � ����� ��������� ��������� �������
	INTEGER N,M,EL_NUMBER,ierr
	!������ � ������� ������������ ������� ����, � �����
	!� ���������
	REAL(KIND=8) LTR(EL_NUMBER),UTR(EL_NUMBER),DIAG(EL_NUMBER),F(N),
	8X(N)	
	!�������, �������� ������� ������� ����
	INTEGER IPTR(N+1),JPTR(EL_NUMBER)

      real*8, allocatable :: W(:)		!������ ��� ��������� ����
      real*8, allocatable :: R0(:)		!������ �������
      real*8, allocatable :: RES1(:)	    !��������������� ������� ��� �������
      real*8, allocatable :: RES2(:)	    !������, ��������� � �������������������
     	real*8, allocatable :: L(:)		!������ � ������� ������������ ������� ����, �
     	real*8, allocatable :: U(:)		!����� � ���������, � ������� ��������� ��������
     	real*8, allocatable :: D(:)		!������������ ���������� (ILU-������������)
	!������� �������� �������	
	real*8, allocatable :: GIVENS_MATRIX(:,:)
	real*8, allocatable :: G(:)!������ ������ �����, ���������� � ������
							!������ GMRES
	real*8, allocatable :: Y(:)	!������������ �������� ����������
							!�������� ��������������� �������,
							!������������ �������� � �������
	!�������� �������, ����� ������� � ��������� ������������
	!����� ����������������� �������� ������ 
	!��������������� ������� � ������� ������������� H,
	!���������� �� ��������� ��������
	real*8, allocatable :: V(:,:)
	real*8, allocatable :: H(:,:)

	REAL(KIND=8) EPS, BETTA, F_NORM, NORM, SCAL_PRODUCT
	!���������� ��� ����������� ������� ������ ���������
	INTEGER H1,M1,S1,SS1,H2,M2,S2,SS2,TM
	REAL(KIND=8) TEMP1,TEMP2
	INTEGER ITER_NUM,I,K,T,P,J, INTIK

	allocate(W(N),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in W'

	allocate(R0(N),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in R0'

	allocate(RES1(N),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in RES1'

	allocate(RES2(N),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in RES2'

	allocate(L(EL_NUMBER),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in L'

	allocate(D(N),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in D'

	allocate(U(EL_NUMBER),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in U'

	allocate(GIVENS_MATRIX(2,M),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in GIVENS_MATRIX'

	allocate(G(M+1),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in G'

	allocate(V(N,M+1),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in V(N,M+1)'

	allocate(H(M+1,M),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in H(M+1,M)'

	allocate(Y(M),stat=ierr)
	if(ierr.ne.0)stop 'Allocating error in Y'

	DO I=1,N
	X(I)=0d0; W(I)=0d0; R0(I)=0d0
	END DO
	Y(1:M)=0d0
	G(1:M+1)=0d0
	GIVENS_MATRIX(1:2,1:M)=0d0
	V(1:N,1:M+1)=0d0
	H(1:M+1,1:M)=0d0
	L(1:EL_NUMBER)=0d0
	U(1:EL_NUMBER)=0d0
	D(1:N)=0d0
	ITER_NUM=0
	INTIK=0
	
	!������ �������� �������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	EPS=1D-6              
	
	!��������� ILU-������������ � ������� �������� �������
	CALL L_Usq(IPTR,JPTR,LTR,UTR,DIAG,L,U,D,INTIK,N,M,EL_NUMBER)


	!�������, ��� ��������� ����������� �����
	!�������� ������� ��������. �����, �� 
	!���������������� ����, ������� ��������
	!� ������ ������ ����������������� �������
	CALL LOW_TRIANGLE_TASK(IPTR,JPTR,D,L,F,R0,N,M,EL_NUMBER)

	!��������� �������� ����� ������� ��� ���������� �����������
	BETTA=NORM(R0,N,M,EL_NUMBER)
	!��������� �������� ����� ��� ������ �����
	F_NORM=NORM(F,N,M,EL_NUMBER)

	DO WHILE(BETTA/F_NORM .GT. EPS)
	
	  !���������� ������ G
	  G(1)=BETTA; G(2:M+1)=0d0
	  !���������� ������, ����������� ������ ��������������� �������
	  DO K=1,N
		V(K,1)=R0(K)/BETTA
	  END DO

	  P=M

	  DO I=1,M

		!����� ������ W ��� ����������������� �������
	CALL UP_TRIANGLE_TASK(IPTR,JPTR,D,U,V(1,I),RES1,N,M,EL_NUMBER)
	CALL MATRIX_ON_VECTOR_MULT(IPTR,JPTR,LTR,UTR,DIAG,RES1,RES2,N,M,
	8EL_NUMBER)
	CALL LOW_TRIANGLE_TASK(IPTR,JPTR,D,L,RES2,W,N,M,EL_NUMBER)

*		CALL MATRIX_ON_VECTOR_MULT(IPTR,JPTR,LTR,UTR,DIAG,V(1,I),W)

		!���������� ������������ �������� �
		!��������� �������� ������
		DO K=1,I
			H(K,I)=SCAL_PRODUCT(W,V(1,K),N,M,EL_NUMBER)
			DO T=1,N
				W(T)=W(T)-H(K,I)*V(T,K)
			END DO	
		END DO
		
		!�������������� ���������� ������
		H(I+1,I)=NORM(W,N,M,EL_NUMBER)
		DO T=1,N
			V(T,I+1)=W(T)/H(I+1,I)
		END DO
		
		!�������� ��� ���������� ������� ��������
		!�� ��������� i-�� ������� ������� H
		DO K=1,I-1
		TEMP1=GIVENS_MATRIX(1,K)*H(K,I) - GIVENS_MATRIX(2,K)*H(K+1,I)
		TEMP2=GIVENS_MATRIX(2,K)*H(K,I) + GIVENS_MATRIX(1,K)*H(K+1,I)
		H(K,I)=TEMP1
		H(K+1,I)=TEMP2
		END DO

		!�������� ������� �������, ������� ��������
		!������� H(i+1,i)
		CALL FORM_GIVENS_ROTATION_MATRIX(H(I,I),H(I+1,I),
     x			  GIVENS_MATRIX(1,I),GIVENS_MATRIX(2,I))
		
		!����� ���������� �� ��������� ������� ��������,
		!������� ������������ ������� ������� H
		H(I,I)=SQRT(H(I,I)*H(I,I) + H(I+1,I)*H(I+1,I))

		!����� ���������� �������� ��������
		!�������� ������ ������� G
		TEMP1=GIVENS_MATRIX(1,I)*G(I) !-GIVENS_MATRIX(2,I)*G(I+1)
		TEMP2=GIVENS_MATRIX(2,I)*G(I) !+GIVENS_MATRIX(1,I)*G(I+1)
		G(I)=TEMP1
		G(I+1)=TEMP2

		!���� ��������� ���������� ������� G ��������� ������,
		!��� �������� ��������, �� ������ ��� ������
		!���������� ����������� �������������� ��������
		!��������������� �������, � ����� ����� ����������
		!� ������� �������, ��� ����������� ��������� ������,
		!��� �������������� ����������
		IF( ABS(G(I+1)) .LT. EPS) THEN
			P=I		!���������� ����� ����������� �������
			EXIT	!������� �� ����� �� I
		END IF
		
	  END DO !(DO I=1,M)

	  !������� ����������������� ���� (PxP) ���
	  !������ ��������� ���� �� �������
	  DO I=P,1,-1
		Y(I)=G(I)
		DO J=I+1,P
			Y(I)=Y(I)-H(I,J)*Y(J)
		END DO
		Y(I)=Y(I)/H(I,I)
	  END DO

	  !� ������������ � ���������� ��������������,
	  !������� ������� �������
	  DO I=1,P
		DO J=1,N
			X(J)=X(J)+Y(I)*V(J,I)
		END DO
	  END DO

	  ITER_NUM=ITER_NUM+1

	  !������� �������, �������������
	  !����� �������� X
	  CALL UP_TRIANGLE_TASK(IPTR,JPTR,D,U,X,RES1,N,M,EL_NUMBER)
	CALL MATRIX_ON_VECTOR_MULT(IPTR,JPTR,LTR,UTR,DIAG,RES1,RES2,N,
	8M,EL_NUMBER)
	  DO I=1,N
		RES2(I)=F(I)-RES2(I)
	  END DO
	  CALL LOW_TRIANGLE_TASK(IPTR,JPTR,D,L,RES2,R0,N,M,EL_NUMBER)

	  BETTA=NORM(R0,N,M,EL_NUMBER)

	  PRINT*,'BETTA=',BETTA/F_NORM,'  P=',P

	  !���� ����������� ���������� ����� ��������
	  !�������� ��������������� ������� �������,
	  !�� ���������� ������ ������
	  IF(P .LT. M) EXIT

	END DO !DO WHILE(BETTA .GT. EPS)
	
	CALL UP_TRIANGLE_TASK(IPTR,JPTR,D,U,X,RES1,N,M,EL_NUMBER)


	!����� ����������� �� �����
	PRINT*,'---------------------------------------------'
	OPEN(UNIT=25,FILE='ANALITIC SOLUTION.TXT')
	OPEN(UNIT=30,FILE='out')!RESULT.TXT
	OPEN(UNIT=50,FILE='GMRES_ERROR.TXT')
	OPEN(UNIT=55,FILE='LOS_ERROR.TXT')
	PRINT*,'DIMMENSION OF PROBLEM=',N
	PRINT*,'DIMMENSION OF KRYLOVS SUBSPACE=',M
	PRINT*,'EPSILON=',EPS
	PRINT*,'FINAL RESIDUAL=',BETTA/F_NORM
	PRINT*,'NUMBER OF ITERATIONS=',ITER_NUM
	PRINT*,'---------------------------------------------'
	PRINT*,'    GMRES ERROR:        ','    LOS ERROR:'
	DO I=1,N
!		READ(25,*) TEMP1
!		READ(30,*) TEMP2
!		PRINT*,TEMP1-X(I), TEMP1-TEMP2
!		PRINT*,X(I)
		WRITE(30,*) RES1(I)
	END DO
	PRINT*,'---------------------------------------------'
	CLOSE(25)
	CLOSE(30)
	CLOSE(50)
	CLOSE(55)
	deallocate(GIVENS_MATRIX)
	deallocate(G)
	deallocate(Y)
	deallocate(V)
	deallocate(H)
	deallocate(W)
	deallocate(R0)
	deallocate(RES1)
	deallocate(RES2)
	deallocate(L)
	deallocate(D)
	END

*****************************************************
	!����� �������
	FUNCTION NORM(X,N,M,EL_NUMBER)
	IMPLICIT NONE
	
	INTEGER N,M,EL_NUMBER,I
	REAL(KIND=8) NORM,TEMP,X(N)
	
	TEMP=0d0
	DO I=1,N
		TEMP=TEMP + X(I)*X(I)
	END DO
	NORM=SQRT(TEMP)

	END

************************************************
	!��������� ������������ ��������
	FUNCTION SCAL_PRODUCT(X,Y,N,M,EL_NUMBER)
	IMPLICIT NONE

	INTEGER N,M,EL_NUMBER,I
	REAL(KIND=8) SCAL_PRODUCT,TEMP,X(N),Y(N)

	TEMP=0d0
	DO I=1,N
		TEMP=TEMP + X(I)*Y(I)
	END DO
	SCAL_PRODUCT=TEMP

	END

************************************************
*	��������� �������, �������� � ����������-��������
*	�������, �� ������ Y
	SUBROUTINE MATRIX_ON_VECTOR_MULT(IPTR,JPTR,LTR,UTR,DIAG,Y,RES,N,M,
	8EL_NUMBER)
	IMPLICIT NONE
	
	INTEGER N,M,EL_NUMBER
	REAL(KIND=8) LTR(*),UTR(*),DIAG(*),Y(N),RES(N)
	INTEGER IPTR(*),JPTR(*)
	INTEGER I,J,K

	DO I=1,N
		RES(I)=DIAG(I)*Y(I)
	END DO

	DO I=1,N
		DO J=IPTR(I),IPTR(I+1)-1
			K=JPTR(J)
			RES(I)=RES(I)+LTR(J)*Y(K)
			RES(K)=RES(K)+UTR(J)*Y(I)

		END DO

	END DO
	END

************************************************
*	���������� ������� �������� �������
	SUBROUTINE FORM_GIVENS_ROTATION_MATRIX(H1,H2,C,S)
	IMPLICIT NONE
	REAL(KIND=8) H1,H2,	!�������� ������������� �������
     x			 C,S,	!�������� ������� �������
     x			 TEMP	!"TEMP" - �� � � ������ "TEMP"
	
	TEMP=SQRT(H1*H1 + H2*H2)
	
	C=H1/TEMP
	S=-H2/TEMP
	
	END

************************************************
*	�������� LU(sq)-������������ ����������
	subroutine L_Usq(ig,jg,ggl,ggu,di,ult,ut,ulu,intik,M_DIM,M1,
	8EL_NUMBER) 
	IMPLICIT NONE
	
	INTEGER M_DIM,M1,EL_NUMBER

      REAL(KIND=8) ggl(*),ggu(*),di(*),ult(*),ut(*),ulu(*)
	REAL(KIND=8) z,z1,z2,k,k1,m,i1
	integer*4 i,j,l,ll,intik,ig(*),jg(*)

	if(intik.eq.1) then
	ulu(1)=sqrt(di(1))
	do i=2,M_DIM
		z2=di(i)
		m=ig(i+1)-ig(i)
		i1=ig(i)
	    do l=i1,i1+m-1
			z=ggl(l)
			z1=ggu(l)
			k=ig(jg(l))     !� ������� �������� � ������
			do ll=i1,l-1
				k1=jg(ll)
				j=0
				do while(jg(k)<=jg(ll).and.j==0)
					if(jg(k)==jg(ll))then
										j=1
					end if
					k=k+1
				end do
				if(j==1)then
					k=k-1
					z=z-ult(ll)*ut(k)
					z1=z1-ut(ll)*ult(k)
				end if
			end do
			ult(l)=z/ulu(jg(l))
			ut(l)=z1/ulu(jg(l))
			z2=z2-ult(l)*ut(l)
		end do
		ulu(i)=sqrt(z2)
	end do
	else
	do i=1,M_DIM
	 ulu(i)=1d0 !sqrt(di(i)) !
	end do
	do i=1,EL_NUMBER 
		ult(i)=0d0;
		ut(i)=0d0;
	end do
	end if  
	end

************************************************
*	������� ���������������� ����
	SUBROUTINE LOW_TRIANGLE_TASK(IPTR,JPTR,D,L,F,V,MATRIX_DIM,M,EL_NUMBER)
	IMPLICIT NONE
	
	INTEGER MATRIX_DIM,M,EL_NUMBER
	REAL(KIND=8) L(*),D(*),F(*),V(*)
	INTEGER IPTR(*),JPTR(*)
	INTEGER I,J

	DO I=1,MATRIX_DIM
		V(I)=F(I)
		DO J=IPTR(I),IPTR(I+1)-1
			V(I)=V(I)-L(J)*V(JPTR(J))
		END DO
		V(I)=V(I)/D(I)
	END DO

	END

************************************************
*	������� ����������������� ����
	SUBROUTINE UP_TRIANGLE_TASK(IPTR,JPTR,D,U,F,V,MATRIX_DIM,M,EL_NUMBER)
	IMPLICIT NONE
	
	INTEGER MATRIX_DIM,M,EL_NUMBER
	REAL(KIND=8) U(*),D(*),F(*),V(*)
	INTEGER IPTR(*),JPTR(*)
	INTEGER I,J
	
	DO I=1,MATRIX_DIM
	V(I)=F(I)
	END DO
	V(MATRIX_DIM)=V(MATRIX_DIM)/D(MATRIX_DIM)

	DO I=MATRIX_DIM,2,-1
		
		DO J=IPTR(I),IPTR(I+1)-1
			V(JPTR(J))=V(JPTR(J))-U(J)*V(I)
		END DO
		V(I-1)=V(I-1)/D(I-1)
	END DO

	END
