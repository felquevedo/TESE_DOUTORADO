*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, pVolDer, hrmflg, var3, var4,
     &                   var5, var6, var7)
#include "impcom.inc"
      !**************************************************************************!
      !** SUBROTINA: USERMAT                                                   **!
      !**                                                                      **!
      !** Objetivo: seleciona as subrotinas que atualizam as tensões,          **!
      !**           variáveis de estado e matriz constitutiva                  **!
      !**                                                                      **!
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!         
      !
      ! OBSERVAÇÕES:
      ! ------------
      !      
      !     - Essa rotina funciona para 3D, EPD e AXI
      ! 
      !**************************************************************************! 
      ! Informações sobre as variáveis                                           !
      !**************************************************************************! 
      ! 
      ! CÓDIGOS:
      ! --------
      ! 
      !     int - inteiro
      !     dp  - dupla precisão
      !     sc  - escalar
      !     ar  - arranjo(número de componentes)
      !     i   - entrada
      !     o   - saída
      !     io  - entrada e saída
      !     l   - local
      ! 
      ! VARIÁVEIS DE ENTRADA:
      ! ----------------------
      ! 
      !  matId      (int,sc,i)               número do material
      !  elemId     (int,sc,i)               número do elemento
      !  kDomIntPt  (int,sc,i)               k-ésimo ponto de integração
      !  kLayer     (int,sc,i)               k-ésima camada
      !  kSectPt    (int,sc,i)               k-ésimo ponto de integração da seção
      !  ldstep     (int,sc,i)               número do passo
      !  isubst     (int,sc,i)               número do subpasso
      !  nDirect    (int,sc,in)              quantidade de componentes normais
      !  nShear     (int,sc,in)              quantidade de componentes de cisalhamento
      !  ncomp      (int,sc,in)              total de componentes (nDirect + nShear)
      !  nstatev    (int,sc,i)               quantidade de variáveis de estado
      !  nProp      (int,sc,i)               quantidade de constantes do material
      !  Temp       (dp,sc,in)               temperatura no inicio do incremento de tempo
      !  dTemp      (dp,sc,in)               incremento de temperatura
      !  Time       (dp,sc,in)               tempo no inicio do incremento (t)
      !  dTime      (dp,sc,in)               incremento de tempo (dt)
      !  Strain     (dp,ar(ncomp),i)         deformação no inicio do incremento de tempo
      !  dStrain    (dp,ar(ncomp),i)         incremento de deformação
      !  prop       (dp,ar(nprop),i)         constantes do material definidas pelo comando TB,USER
      !  coords     (dp,ar(3),i)             coordenadas
      !  defGrad_t  (dp,ar(3,3),i)           gradiente de deformação no tempo t
      !  defGrad    (dp,ar(3,3),i)           gradiente de deformação no tempo t+dt
      !  hrmflg     (dp,sc,io)               flag que indica análise harmônica (não usada)
      ! 
      ! VARIÁVEIS DE ENTRADA E SAIDA:
      ! ------------------------------
      !        
      !  stress     (dp,ar(ncomp),io)        tensão
      !  ustatev    (dp,ar(nstatev),io)      variaveis de estado do usuario
      !  sedEl      (dp,sc,io)               trabalho elástico
      !  sedPl      (dp,sc,io)               trabalho plástico
      !  epseq      (dp,sc,io)               deformação plástica equivalente
      !  epsPl      (dp,ar(ncomp),io)        deformação plastica
      !  var?       (dp,sc,io)               variáveis não utilizadas
      ! 
      ! VARIÁVEIS DE SAIDA:
      ! --------------------
      ! 
      !  keycut     (int,sc,io)             inteiro que controla o método bisect/cut (não usada)
      !                                     (variável que será determinada no controle de solução do ANSYS)
      !  dsdePl     (dp,ar(ncomp,ncomp),io) matriz constitutiva tangente (ou algoritmica)
      !  pVolDer    (dp,ar(3),o)            derivada do potencial volumétrico com relação a J (não usada)
      !                                         pVolDer(1) = dU/dJ
      !                                         pVolDer(2) = d^2U/dJ^2
      !                                         pVolDer(3) = d^3U/dJ^3
      !  tsstif     (dp,ar(2),o)            módulo de rigidez transversal
      !                                         tsstif(1) - Gxz
      !                                         tsstif(2) - Gyz
      !                                         tsstif(1) é também utilizado para calcular rigidez hourglass 
      !                                         este valor deve ser definido quando elementos de baixa ordem
      !                                         tais como 181, 182 e 185 com integração uniforme for usada.
      !  epsZZ      (dp,sc,o)               deformação ZZ para estado plano de tensões (não usada)
      !  cutFactor  (dp,sc,o)               fator de corte do tamanho do passo (recomendado entre 0 e 1) (não usada)  
      ! 
      ! VARIÁVEIS LOCAIS:
      ! -----------------
      !  
      !  modelo     (int,sc,l)              identifica o modelo do material
      !
      !**************************************************************************! 
      ! Informações sobre as componentes e notação de Voigt adotada              !
      !**************************************************************************! 
      ! 
      !  ncomp   6   para 3D  (nshear=3)
      !  ncomp   4   para EPD e AXI (nShear = 1)
      ! 
      !  Vetores de tensão e deformação
      ! 
      !      11, 22, 33, 12, 23, 13    para 3D
      !      11, 22, 33, 12            para EPD e AXI
      ! 
      !  Matriz constitutiva 3D
      !       dsdePl    |  1111   1122   1133   1112   1123   1113 |
      !       dsdePl    |  2211   2222   2233   2212   2223   2213 |
      !       dsdePl    |  3311   3322   3333   3312   3323   3313 |
      !       dsdePl    |  1211   1222   1233   1212   1223   1213 |
      !       dsdePl    |  2311   2322   2333   2312   2323   2313 |
      !       dsdePl    |  1311   1322   1333   1312   1323   1313 |
      !  
      !   Matriz constitutiva EPD e AXI (11, 22, 33, 12)
      !       dsdePl    |  1111   1122   1133   1112 |
      !       dsdePl    |  2211   2222   2233   2212 |
      !       dsdePl    |  3311   3322   3333   3312 |
      !       dsdePl    |  1211   1222   1233   1212 |
      ! 
      !**************************************************************************! 
      ! Declaração das variáveis                                                 !
      !**************************************************************************!  
      !
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp,
     &                 modmat
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 pVolDer (3),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION hrmflg
      !
      !
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      DATA             var1/0.0d0/
      DATA             var2/0.0d0/
      ! 
      !**************************************************************************! 
      ! Selecionando subrotinas do material                                      !
      !**************************************************************************! 
      !
      ! modmat = 4 - elastoplástico von-Mises bilinear (nativa do ANSYS)
      ! modmat = 5 - usermat3D_elastico
      ! modmat = 6 - usermat3D_EP
      !
      modmat = prop(1)
      !
      SELECT CASE (modmat)
      !
      CASE(4) 
              CALL usermat3d_vm_ansys (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7) 
      CASE(5)
              CALL usermat3d_elastico (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)   
      CASE(6)
              CALL usermat3d_EP (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)                            
      END SELECT     
      !
      ! 
      !
      !
      RETURN
      END
      ! 
      !
      !
      ! --------------------------------------------------------------------------
      !
*deck,usermat3d    USERDISTRIB  parallel                                gal
      SUBROUTINE usermat3d_elastico (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
#include "impcom.inc"
      !**************************************************************************!
      !** SUBROTINA: USERMAT3D_ELASTICO                                        **!
      !**                                                                      **!
      !** Objetivo: atualiza as tensões, variáveis de estado e matriz          **!
      !**           constitutiva para o modelo constitutivo elástico           **!
      !**                                                                      **!
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!         
      ! 
      !**************************************************************************! 
      ! Declaração variáveis de entrada e saída da subrotina                     !
      !**************************************************************************! 
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      ! 
      !**************************************************************************! 
      ! Informações sobre as variáveis locais                                    !
      !**************************************************************************! 
      ! i,j           (int,sc,l)               contadores
      ! mcomp         (int,sc,l)               total de componentes
      ! stresstrial   (dp,ar(6),l)             tensão tentativa
      ! dsdeEl        (dp,ar(6,6),l)           matriz constitutiva elastica
      ! young         (dp,sc,l)                módulo de Young
      ! posn          (dp,sc,l)                coeficiente de Poisson
      ! 
      !**************************************************************************! 
      ! Declaração das variáveis locais                                          !
      !**************************************************************************! 
      INTEGER             i, j, mcomp
      PARAMETER           (mcomp = 6)
      DOUBLE PRECISION    dstress0(mcomp), dsdeEl(mcomp,mcomp),
     &                    sigi(mcomp), young, posn, stresstrial(mcomp)
      EXTERNAL            vzero, vmove, get_ElmData, get_ElmInfo
      ! 
      !**************************************************************************! 
      ! Entrada de dados                                                         !
      !**************************************************************************! 
      young    = prop(2)
      posn     = prop(3)
      ! 
      !**************************************************************************! 
      ! Calculo da matriz constitutiva                                           !
      !**************************************************************************! 
      dsdeEl   = 0.0d0
      CALL MatrizD(young,posn,dsdeEl)
      ! 
      !**************************************************************************! 
      ! Calculo módulo de rigidez transversal para hourglass                     !
      !**************************************************************************!  
      tsstif(1) = 0.5d0*(young /(1.0d0+posn))
      !
      !**************************************************************************! 
      ! Calculo preditor elástico                                                !
      !**************************************************************************! 
      stresstrial = 0.0d0
      stresstrial = MATMUL(dsdeEl(1:ncomp,1:ncomp),
     &              strain(1:ncomp)+dstrain(1:ncomp))
      !
      !**************************************************************************! 
      ! Atualizando a tensão e variáveis de estado                               !
      !**************************************************************************!       
      stress = stresstrial
      !
      !**************************************************************************! 
      ! Atualizando o módulo constitutivo                                        !
      !**************************************************************************!       
      dsdePl = dsdeEl
      !
      RETURN
      END
      !
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d_EP (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
#include "impcom.inc"  
      !**************************************************************************!
      !** SUBROTINA: USERMAT3D_EP                                              **!
      !**                                                                      **!
      !** Objetivo: atualiza as tensões, variáveis de estado e matriz          **!
      !**           constitutiva para o modelo constitutivo elastoplastico     **!
      !**                                                                      **!
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!         
      ! 
      !**************************************************************************! 
      ! Declaração variáveis de entrada e saída da subrotina                     !
      !**************************************************************************! 
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      ! 
      !**************************************************************************! 
      ! Informações sobre as variáveis locais                                    !
      !**************************************************************************! 
      ! i,j           (int,sc,l)               contadores
      ! mcomp         (int,sc,l)               total de componentes
      ! stresstrial   (dp,ar(6),l)             tensão tentativa
      ! dsdeEl        (dp,ar(6,6),l)           matriz constitutiva elastica
      ! young         (dp,sc,l)                módulo de Young
      ! posn          (dp,sc,l)                coeficiente de Poisson
      ! 
      !**************************************************************************! 
      ! Declaração das variáveis locais                                          !
      !**************************************************************************! 
      INTEGER             i, j                ! contadores
      INTEGER             mcomp               ! máximo de componentes
      INTEGER             nrmax               ! máximo número de interações de NR
      DOUBLE PRECISION    vPi                 ! Pi
      PARAMETER           (mcomp = 6)         ! máximo de componentes
      PARAMETER           (nrmax = 10)        ! máximo número de interações de NR
      PARAMETER           (vPi = 3.14159265358979323846d0) ! máximo número de interações de NR
      DOUBLE PRECISION    dstress(mcomp)      ! incremento de tensão
      DOUBLE PRECISION    dsdeEl(mcomp,mcomp) ! matriz constitutiva elástica
      DOUBLE PRECISION    dsdeElinv(mcomp,mcomp) ! inversa da matriz constitutiva elástica
      DOUBLE PRECISION    sigi(mcomp)         ! tensões iniciais
      !
      DOUBLE PRECISION    young               ! modulo de Young
      DOUBLE PRECISION    posn                ! coeficiente de Poisson
      INTEGER             superficief         ! funçaõ de escoamento: 1- (MC ou TR), 2- (DP ou VM)
      INTEGER             superficieg         ! função potencial: 1- (MC ou TR), 2- (DP ou VM)
      INTEGER             flag                ! 1- DP coincide com as bordas externas MC 
      DOUBLE PRECISION    c                   ! coesão
      DOUBLE PRECISION    fi                  ! angulo de atrito
      DOUBLE PRECISION    psi                 ! angulo de dilatância
      !
      DOUBLE PRECISION    stresstrial(mcomp)  ! tensão tentativa
      DOUBLE PRECISION    I1                  ! primeiro invariante das tensões
      DOUBLE PRECISION    J2                  ! segundo invariante do desviador
      DOUBLE PRECISION    J3                  ! terceiro invariante do desviador
      DOUBLE PRECISION    p                   ! pressão hidrostática
      DOUBLE PRECISION    q                   ! tensão equivalente
      DOUBLE PRECISION    theta               ! angulo de Lode
      DOUBLE PRECISION    s(mcomp)            ! tensor desviador
      DOUBLE PRECISION    g1(mcomp),g2(mcomp),g3(mcomp)   ! direções do vetor de fluxo
      DOUBLE PRECISION    c1,c2,c3            ! magnitudes do vetor de fluxo
      DOUBLE PRECISION    dgds(mcomp)         ! vetor de fluxo
      DOUBLE PRECISION    dfds(mcomp)         ! derivada de f em relação as tensões
      DOUBLE PRECISION    f                   ! função de escoamento
      DOUBLE PRECISION    ftDg                ! dfds*D*dgds
      DOUBLE PRECISION    dlam                ! incremento da magnitude da deformação plástica
      DOUBLE PRECISION    dlamk               ! incremento dlam pelo NR do corretor plástico
      DOUBLE PRECISION    DgftD(mcomp,mcomp)  ! Produto para calcular matriz algoritmica
      DOUBLE PRECISION    dfds_m(1,mcomp)     ! Produto para calcular matriz algoritmica
      DOUBLE PRECISION    dgds_m(mcomp,1)     ! Produto para calcular matriz algoritmica
      
      EXTERNAL            vzero, vmove, get_ElmData, get_ElmInfo
      ! 
      !**************************************************************************! 
      ! Entrada de dados                                                         !
      !**************************************************************************! 
      keycut      = 0
      young       = prop(2)
      posn        = prop(3)
      superficief = prop(4)
      superficieg = prop(5)
      flag        = prop(6)
      c           = prop(7)
      fi          = prop(8)*vPi/180
      psi         = prop(9)*vPi/180
      ! 
      !**************************************************************************! 
      ! Calculo da matriz constitutiva                                           !
      !**************************************************************************! 
      dsdeEl   = 0.0d0
      CALL MatrizD(young,posn,dsdeEl)
      ! 
      !**************************************************************************! 
      ! Calculo módulo de rigidez transversal para hourglass                     !
      !**************************************************************************!  
      tsstif(1) = 0.5d0*(young /(1.0d0+posn))
      !
      !**************************************************************************! 
      ! Calculo preditor elástico                                                !
      !**************************************************************************! 
      stresstrial = 0.0d0
      stresstrial = MATMUL(dsdeEl(1:ncomp,1:ncomp),
     &              strain(1:ncomp)+dstrain(1:ncomp))
      !
      stress = stresstrial
      dsdePl = dsdeEl
      !
      !**************************************************************************! 
      ! Calcula dgds                                                             !
      !**************************************************************************!  
      ! Calcula os invariantes do tensor de tensões
      CALL invars(stress,I1,J2,J3,p,q,theta,s)    
      !
      ! Calcula as direções do vetor de fluxo
      CALL g1g2g3(s,J2,g1,g2,g3)
      !
      ! Calcula as magnitudes componentes do vetor de fluxo
      CALL c1c2c3(J2,theta,superficieg,psi,flag,c1,c2,c3)
      !
      ! Calcula dgds
      dgds = c1*g1 + c2*g2 + c3*g3
      !
      !**************************************************************************! 
      ! Calcula função de escoamento                                             !
      !**************************************************************************!  
      CALL yield(superficief,I1,J2,theta,c,fi,flag,f)
      !
      !**************************************************************************! 
      ! Verifica o critério de escoamento                                        !
      !**************************************************************************! 
      IF(f.GT.0.0d0)THEN
          ! Aplica o corretor plástico
          dlam = 0.0d0
          i = 0
          DO 
              ! Calcula dfds
              CALL invars(stress,I1,J2,J3,p,q,theta,s)
              CALL g1g2g3(s,J2,g1,g2,g3)
              CALL c1c2c3(J2,theta,superficief,fi,flag,c1,c2,c3)
              dfds    = c1*g1 + c2*g2 + c3*g3
              !
              ! Calcula dlamda
              ftDg    = DOT_PRODUCT(dfds,MATMUL(dsdeEl,dgds))
              dlamk   = f/(ftDg)
              dstress = -dlamk*MATMUL(dsdeEl,dgds)
              !
              ! Atualizando variáveis internas
              stress  = stress + dstress
              dlam    = dlam + dlamk
              !
              ! Testa o critério de escoamento
              CALL invars(stress,I1,J2,J3,p,q,theta,s) 
              CALL yield(superficief,I1,J2,theta,c,fi,flag,f)   
              IF(f.LE.0.00001d0)EXIT
              !
              ! Caso atinja o número de iterações limites, faça a bisseção
              i = i + 1
              IF(i.EQ.nrmax)THEN
                  keycut = 1
                  RETURN
              ENDIF
          ENDDO
          !
          ! Calcula as deformações plásticas 
          epsPl = dlam*dgds    
          !
          !**************************************************************************! 
          ! Atualizando o módulo constitutivo                                        !
          !**************************************************************************!       
          dfds_m(1,:) = dfds
          dgds_m(:,1) = dgds
          DgftD = MATMUL(MATMUL(MATMUL(dsdeEl,dgds_m),dfds_m),dsdeEl)
          dsdePl = dsdeEl - DgftD/(ftDg)
          !
      ENDIF
      !
      !
      !
      RETURN
      END
      !
      !
      !
      SUBROUTINE matrizD(E,Poisson,D)
      !**************************************************************************!
      !** Função: matrizD                                                      **!
      !**                                                                      **!
      !** Objetivo: calcula a matriz consitutiva do material isotrópico        **!
      !**                                                                      **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      IMPLICIT NONE
      DOUBLE PRECISION E          ! módulo de elasticidade                  
      DOUBLE PRECISION Poisson    ! coeficiente de Poisson
      DOUBLE PRECISION D(6,6)     ! matriz constitutiva elástica isotrópica
      !
	D=0.0d0
	D(1,1)=(E*(1.0d0-Poisson))/((1.0d0+Poisson)*(1.0d0-2.0d0*Poisson))
	D(1,2)=(E*Poisson)/((1.0d0+Poisson)*(1.0d0-2.0d0*Poisson))
	D(1,3)=D(1,2)
	D(2,1)=D(1,2)
	D(2,2)=D(1,1)
	D(2,3)=D(1,2)
	D(3,1)=D(1,3)
	D(3,2)=D(2,3)
	D(3,3)=D(1,1)
	D(4,4)=(E)/((1.0d0+Poisson)*2.0d0)
	D(5,5)=D(4,4)
	D(6,6)=D(4,4)
      !
      END SUBROUTINE MatrizD
      !
      SUBROUTINE invars(stress,I1,J2,J3,p,q,theta,s)      
      !**************************************************************************!
      !** Subrotina: invars                                                    **!
      !**                                                                      **!
      !** Objetivo: calcula os invariantes do tensor de tensões e o tensor     **!
      !**           desviador                                                  **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!   
      IMPLICIT NONE
      DOUBLE PRECISION stress(6)  ! tensões
      DOUBLE PRECISION I1         ! primeiro invariante do tensor de tensões
      DOUBLE PRECISION J2,J3      ! primeiro e segundo invariante do desviador
      DOUBLE PRECISION s(6)       ! desviador
      DOUBLE PRECISION p,q        ! pressão hidrostática e tensão eq. de vm
      DOUBLE PRECISION theta      ! angulo de Lode
      DOUBLE PRECISION sine       ! variavel auxiliar
      !
      ! Inicializando variaveis
      I1      = 0.0d0
      J2      = 0.0d0
      p       = 0.0d0
      q       = 0.0d0
      s       = 0.0d0
      J3      = 0.0d0
      theta   = 0.0d0
      sine    = 0.0d0
      !
      ! Calculo do I1
      I1 = stress(1) + stress(2) + stress(3)
      !
      ! Calculo do J2
      J2 = 1/6.0d0*((stress(1)-stress(2))**2+(stress(2)-stress(3))**2+
     &    (stress(3)-stress(1))**2)+
     &    stress(4)**2+stress(5)**2+stress(6)**2
      !
      ! Calculo do p
      p = 1/3.0d0*I1
      !
      ! Calculo do desviador s
      s(1) = stress(1) - p
      s(2) = stress(2) - p
      s(3) = stress(3) - p
      s(4) = stress(4)
      s(5) = stress(5)
      s(6) = stress(6)
      !
      ! Calculo do J3
      J3 = s(1)*s(2)*s(3)-s(1)*s(5)*s(5)-s(2)*s(6)*s(6)-s(3)*s(4)*s(4)+
     & 2.0d0*s(4)*s(5)*s(6)
      !
      ! Calculo da tensão equivalente
      q = SQRT(3.0d0*J2)
      !
      ! Calculo do ângulo de lode
      IF(q < 1.E-10)THEN
          theta = 0.0d0
      ELSE
          sine = -3.0d0*SQRT(3.0d0)*J3/(2.0d0*SQRT(J2)**3)
          IF(sine>=1.0d0)sine=1.0d0
          IF(sine<-1.0d0)sine=-1.0d0
          theta=ASIN(sine)/3.0d0
      END IF
      !
      END SUBROUTINE
      !
      SUBROUTINE g1g2g3(s,J2,g1,g2,g3)
      !**************************************************************************!
      !** Subrotina: g1g2g3                                                    **!
      !**                                                                      **!
      !** Objetivo: calcula as direções do vetor de fluxo                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      IMPLICIT NONE
      DOUBLE PRECISION s(6)                   ! desviador
      DOUBLE PRECISION J2                     ! segundo invariante do desviador
      DOUBLE PRECISION g1(6), g2(6), g3(6)    ! direções do vetor de fluxo
      !
      ! Inicializando variaveis
      g1 = 0.0d0
      g2 = 0.0d0
      g3 = 0.0d0
      !
      ! Calculo do g1
      g1(1) = 1.0d0
      g1(2) = 1.0d0
      g1(3) = 1.0d0
      g1(4) = 0.0d0
      g1(5) = 0.0d0
      g1(6) = 0.0d0
      !
      ! Calculo do g2
      g2(1) = s(1)
      g2(2) = s(2)
      g2(3) = s(3)
      g2(4) = 2.0d0*s(4)
      g2(5) = 2.0d0*s(5)
      g2(6) = 2.0d0*s(6)
      g2 = 1/(2*SQRT(J2))*g2
      !
      ! Calculo do g3
      g3(1) = s(2)*s(3) - s(5)**2 + J2/3.0d0
      g3(2) = s(1)*s(3) - s(6)**2 + J2/3.0d0
      g3(3) = s(1)*s(2) - s(4)**2 + J2/3.0d0
      g3(4) = 2.0d0*(s(5)*s(6)-s(3)*s(4))
      g3(5) = 2.0d0*(s(6)*s(4)-s(1)*s(5))
      g3(6) = 2.0d0*(s(4)*s(5)-s(2)*s(6))     
      !
      END SUBROUTINE
      !
      SUBROUTINE c1c2c3(J2,theta,superficie,fi,flag,c1,c2,c3)
      !**************************************************************************!
      !** Subrotina: c1c2c3                                                    **!
      !**                                                                      **!
      !** Objetivo: calcula a magnitude das componentes do vetor de fluxo      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      IMPLICIT NONE
      DOUBLE PRECISION J2             ! segundo invariante do desviador
      DOUBLE PRECISION theta          ! ângulo de Lode
      DOUBLE PRECISION fi             ! Ângulo de atrito
      INTEGER          superficie     ! 1- (MC ou TR), 2- (DP ou VM)
      INTEGER          flag           ! 1- DP coincide com as bordas externas MC
      DOUBLE PRECISION c1,c2,c3       ! magnitude das componentes do vetor
      DOUBLE PRECISION cos3th         ! variavel auxiliar
      DOUBLE PRECISION beta           ! parâmetros do DP
      !
      ! Seleciona o modelo
      SELECT CASE (superficie)
      CASE(1)
          !
          ! Mohr-Coulomb ou Tresca (fi = 0)
          !
          c1 = 1.0d0/3.0d0*SIN(fi)
          !
          ! Teste da singularidade
          cos3th = COS(3.0d0*theta)
          IF (cos3th < 1.E-10)THEN
              c2 = 1.0d0
              c3 = 0.0d0
          ELSE   
              c2 = COS(theta)*(1.0d0+TAN(theta)*TAN(3.0d0*theta)+
     &            (SIN(fi)*(TAN(3*theta)-TAN(theta)))/(SQRT(3.0d0)))
              c3 = (SQRT(3.0d0)*SIN(theta)+COS(theta)*SIN(theta))/
     &            (2.0d0*J2*COS(3*theta))
          ENDIF
          CASE(2)
          !
          ! Drucker-Prager ou von-Mises (fi = 0)
          !
          IF(flag.EQ.1)THEN
              ! DP coincide com as bordas mais externas de MC
              beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0-SIN(fi)))
          ELSE
              ! DP coincide com as bordas mais internas de MC
              beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0+SIN(fi)))
          ENDIF
          c1 = beta
          c2 = 1.0d0
          c3 = 0.0d0   
      END SELECT
      !
      END SUBROUTINE
      !
      SUBROUTINE yield(superficie,I1,J2,theta,c,fi,flag,f)
      !**************************************************************************!
      !** Subrotina: c1c2c3                                                    **!
      !**                                                                      **!
      !** Objetivo: calcula o critério de escoamento                           **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      IMPLICIT NONE
      INTEGER          superficie     ! 1-(MC ou TR), 2-(DP ou VM)
      DOUBLE PRECISION I1             ! primeiro invariante do tensor de tensões
      DOUBLE PRECISION J2             ! segundo invariante do desviador
      DOUBLE PRECISION theta          ! ângulo de Lode
      DOUBLE PRECISION c              ! coesão
      DOUBLE PRECISION fi             ! Ângulo de atrito
      INTEGER          flag           ! 1- DP coincide com as bordas externas MC
      DOUBLE PRECISION f              ! função de escoamento
      DOUBLE PRECISION beta,k         ! Parametros para DP
      !
      ! Seleciona o modelo
      SELECT CASE (superficie)
      CASE(1)
          !
          ! Mohr-Coulomb ou Tresca (fi = 0)
          !
          f = SIN(fi)/3.0d0*I1 +  
     &        (COS(theta)-1.0D0/SQRT(3.0d0)*SIN(theta)*SIN(fi))*SQRT(J2)
     &        -c*COS(fi)
          !
          CASE(2)
          !
          ! Drucker-Prager ou von-Mises (fi = 0)
          !
          IF(flag.EQ.1)THEN
              ! DP coincide com as bordas mais externas de MC
              beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0-SIN(fi)))
              k    = 6.0d0*COS(fi)/(SQRT(3.0d0)*(3.0d0-SIN(fi)))*c
          ELSE
              ! DP coincide com as bordas mais internas de MC
              beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0+SIN(fi)))
              k    = 6.0d0*COS(fi)/(SQRT(3.0d0)*(3.0d0+SIN(fi)))*c
          ENDIF
          f = beta*I1 + SQRT(J2) - k
      END SELECT
      !
      END SUBROUTINE
      !
      SUBROUTINE invert(A,ndim)      
      !**************************************************************************!
      !** Subrotina: invart                                                    **!
      !**                                                                      **!
      !** Objetivo: inverte uma matriz pequena                                 **!
      !**                                                                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!   
      IMPLICIT NONE
      INTEGER          ndim                   ! numero de dimensões
      DOUBLE PRECISION A(ndim,ndim)           ! matriz a ser invertida
      DOUBLE PRECISION con,i,k              ! variável auxiliar
      !
      DO k=1,ndim
        con=A(k,k)
        A(k,k)=1.0d0
        A(k,:)=A(k,:)/con
        DO i=1,ndim
          IF(i/=k)THEN
            con=A(i,k)
            A(i,k)=1.0d0
            A(i,:)=A(i,:)-A(k,:)*con
          END IF
        END DO
      END DO
      !
      END SUBROUTINE
      !
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d_vm_ansys (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c      
c
c*************************************************************************
c    INFORMAÇÃO VARIÁVEIS LOCAIS
c*************************************************************************
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,l)            three time of shear moduli
c
c --- temperary variables for solution purpose
c      i, j
c      threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
c  
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
c
      EXTERNAL         vzero, vmove, get_ElmData, get_ElmInfo
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j, ncompgt
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3)
c*************************************************************************
c
      keycut   = 0
      cutFactor = 0.d0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(2)
      posn     = prop(3)
      sigy0    = prop(4)
c *** plastic strain tensor
      call vmove(ustatev(2), epsPl(1), ncomp)
c *** calculate plastic slope
      dsigdep  = young*prop(5)/(young-prop(5))
      twoG     = young / (ONE+posn)
      threeG   = ONEHALF * twoG
      elast1=young*posn/((1.0D0+posn)*(1.0D0-TWO*posn))
      elast2=HALF*twoG
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = elast2
c
c *** calculate elastic stiffness matrix (3d)
c
      dsdeEl(1,1)=(elast1+TWO*elast2)*G(1)*G(1)
      dsdeEl(1,2)=elast1*G(1)*G(2)+elast2*TWO*G(4)*G(4)
      dsdeEl(1,3)=elast1*G(1)*G(3)+elast2*TWO*G(5)*G(5)
      dsdeEl(1,4)=elast1*G(1)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(1,5)=elast1*G(1)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(1,6)=elast1*G(1)*G(6)+elast2*TWO*G(4)*G(5)
      dsdeEl(2,2)=(elast1+TWO*elast2)*G(2)*G(2)
      dsdeEl(2,3)=elast1*G(2)*G(3)+elast2*TWO*G(6)*G(6)
      dsdeEl(2,4)=elast1*G(2)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(2,5)=elast1*G(2)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(2,6)=elast1*G(2)*G(6)+elast2*TWO*G(2)*G(6)
      dsdeEl(3,3)=(elast1+TWO*elast2)*G(3)*G(3)
      dsdeEl(3,4)=elast1*G(3)*G(4)+elast2*TWO*G(5)*G(6)
      dsdeEl(3,5)=elast1*G(3)*G(5)+elast2*TWO*G(5)*G(3)
      dsdeEl(3,6)=elast1*G(3)*G(6)+elast2*TWO*G(6)*G(3)
      dsdeEl(4,4)=elast1*G(4)*G(4)+elast2*(G(1)*G(2)+G(4)*G(4))
      dsdeEl(4,5)=elast1*G(4)*G(5)+elast2*(G(1)*G(6)+G(5)*G(4))
      dsdeEl(4,6)=elast1*G(4)*G(6)+elast2*(G(4)*G(6)+G(5)*G(2))
      dsdeEl(5,5)=elast1*G(5)*G(5)+elast2*(G(1)*G(3)+G(5)*G(5))
      dsdeEl(5,6)=elast1*G(5)*G(6)+elast2*(G(4)*G(3)+G(5)*G(6))
      dsdeEl(6,6)=elast1*G(6)*G(6)+elast2*(G(2)*G(3)+G(6)*G(6))
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
c
c
c *** get initial stress
      call get_ElmInfo('NCOMP', ncompgt)
      call vzero(sigi(1),ncompgt)
      call get_ElmData ('ISIG', elemId,kDomIntPt, ncompgt, sigi)
c
c *** calculate the trial stress and
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         strainEl(i) = Strain(i) + dStrain(i)
      end do
      call vzero(sigElp, 6)
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i)+dsdeEl(j,i)*strainEl(j)
         end do
         sigElp(i) = sigElp(i)
      end do


c *** Update stress in case of elastic/unloading
      do i=1,ncomp
         stress(i) = sigElp(i)
      end do


c
!      cutFactor = 0.5d0
      return
      end
