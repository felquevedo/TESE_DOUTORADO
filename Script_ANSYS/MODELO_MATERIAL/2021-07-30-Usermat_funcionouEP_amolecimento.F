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
      ! modmat = 7 - usermat3D_EPVP
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
c      CASE(7)
c              CALL usermat3d_EPVP (
c     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
c     &                   ldstep,isubst,keycut,
c     &                   nDirect,nShear,ncomp,nStatev,nProp,
c     &                   Time,dTime,Temp,dTemp,
c     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
c     &                   Strain,dStrain, epsPl, prop, coords,
c     &                   var0, defGrad_t, defGrad,
c     &                   tsstif, epsZZ, cutFactor, 
c     &                   var1, var2, var3, var4, var5,
c     &                   var6, var7)    
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
      CALL MatrizD(young,posn,ncomp,dsdeEl)
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
      INTEGER             nrmax               ! máximo número de interações de NR
      DOUBLE PRECISION    vPi                 ! Pi
      PARAMETER           (nrmax = 10000)        ! máximo número de interações de NR
      PARAMETER           (vPi = 3.14159265358979323846d0) ! máximo número de interações de NR
      INTEGER             ncompgt             ! numero de componentes para o vetor de tensões iniciais
      DOUBLE PRECISION    dstress(ncomp)      ! incremento de tensão
      DOUBLE PRECISION    dsdeEl(ncomp,ncomp) ! matriz constitutiva elástica
      DOUBLE PRECISION    dsdeElinv(ncomp,ncomp) ! inversa da matriz constitutiva elástica
      DOUBLE PRECISION    sigi(ncomp)         ! tensões iniciais
      !
      DOUBLE PRECISION    young               ! modulo de Young
      DOUBLE PRECISION    posn                ! coeficiente de Poisson
      INTEGER             superficief         ! funçaõ de escoamento: 1- (MC ou TR), 2- (DP ou VM)
      INTEGER             superficieg         ! função potencial: 1- (MC ou TR), 2- (DP ou VM)
      INTEGER             flag                ! 1- DP coincide com as bordas externas MC 
      DOUBLE PRECISION    c                   ! coesão
      DOUBLE PRECISION    fi                  ! angulo de atrito
      DOUBLE PRECISION    psi                 ! angulo de dilatância
      DOUBLE PRECISION    ci,cp,cr            ! coesão inicial, pico e residual
      DOUBLE PRECISION    eps1,eps2,eps3      ! deformações plásticas equivalentes limites
      !
      DOUBLE PRECISION    stresstrial(ncomp)  ! tensão tentativa
      DOUBLE PRECISION    stressi(ncomp)      ! tensão da interação anterior
      DOUBLE PRECISION    I1                  ! primeiro invariante das tensões
      DOUBLE PRECISION    J2                  ! segundo invariante do desviador
      DOUBLE PRECISION    J3                  ! terceiro invariante do desviador
      DOUBLE PRECISION    p                   ! pressão hidrostática
      DOUBLE PRECISION    q                   ! tensão equivalente
      DOUBLE PRECISION    theta               ! angulo de Lode
      DOUBLE PRECISION    s(ncomp)            ! tensor desviador
      DOUBLE PRECISION    g1(ncomp),g2(ncomp),g3(ncomp)   ! direções do vetor de fluxo
      DOUBLE PRECISION    f1(ncomp),f2(ncomp),f3(ncomp)   ! direções de dfds
      DOUBLE PRECISION    c1,c2,c3            ! magnitudes do vetor de fluxo
      DOUBLE PRECISION    dgds(ncomp)         ! vetor de fluxo
      DOUBLE PRECISION    dfds(ncomp)         ! derivada de f em relação as tensões
      DOUBLE PRECISION    f                   ! função de escoamento
      DOUBLE PRECISION    ftDg                ! dfds*D*dgds
      DOUBLE PRECISION    dlam                ! incremento da magnitude da deformação plástica
      DOUBLE PRECISION    ddlam               ! incremento dlam pelo NR do corretor plástico
      DOUBLE PRECISION    epsEP(ncomp)        ! deformações plásticas
      DOUBLE PRECISION    depsEP(ncomp)       ! incremento de deformação plástica
      DOUBLE PRECISION    epsEPeq             ! deformação plástica equivalente
      DOUBLE PRECISION    depsEPeq            ! incremento de deformação plástica equivalente
      DOUBLE PRECISION    czao                ! constante para o calculo da deformação plástica equivalente
      DOUBLE PRECISION    DgftD(ncomp,ncomp)  ! Produto para calcular matriz algoritmica
      DOUBLE PRECISION    dfds_m(1,ncomp)     ! Produto para calcular matriz algoritmica
      DOUBLE PRECISION    dgds_m(ncomp,1)     ! Produto para calcular matriz algoritmica
      DOUBLE PRECISION    dhdq,dfdq,dfdc,dcde ! Produtos referente ao endurecimento/amolecimento
      !
      EXTERNAL            vzero, vmove, get_ElmData, get_ElmInfo
      EXTERNAL            g1g2g3
      ! 
      !**************************************************************************! 
      ! Entrada de dados                                                         !
      !**************************************************************************! 
      !
      ! Variavel que controla o método da bisseção caso não haja convergência
      keycut      = 0
      !
      ! Propriedades elásticas
      young       = prop(2)
      posn        = prop(3)
      !
      ! Propriedades do modelo EP
      superficief = prop(4)
      superficieg = prop(5)
      flag        = prop(6)
      fi          = prop(7)*vPi/180
      psi         = prop(8)*vPi/180
      ci          = prop(9)
      cp          = prop(10)
      cr          = prop(11)
      eps1        = prop(12)
      eps2        = prop(13)
      eps3        = prop(14)
      !
      ! Limpando variáveis
      depsEP      = 0.0d0
      c           = 0.0d0
      f           = 0.0d0
      !
      !**************************************************************************! 
      ! Coletando variáveis de estado e deformações plásticas                    !
      !**************************************************************************! 
      epsEP       = 0.0d0
      epsEPeq     = 0.0d0
      stressi     = 0.0d0
      epsEPeq = ustatev(1)
      CALL vmove(ustatev(2), epsEP(1), ncomp)
      CALL vmove(ustatev(8), stressi(1), ncomp)
      ! 
      !**************************************************************************! 
      ! Calculo da matriz constitutiva                                           !
      !**************************************************************************! 
      dsdeEl      = 0.0d0
      CALL MatrizD(young,posn,ncomp,dsdeEl)
      ! 
      !**************************************************************************! 
      ! Calculo módulo de rigidez transversal para hourglass                     !
      !**************************************************************************!  
      tsstif(1) = 0.5d0*(young /(1.0d0+posn))
      !
      !**************************************************************************! 
      ! Coletando tensões iniciais                                               !
      !**************************************************************************! 
      call get_ElmInfo('NCOMP', ncompgt)
      call vzero(sigi(1),ncompgt)
      call get_ElmData ('ISIG', elemId,kDomIntPt, ncompgt, sigi)  
      !
      !**************************************************************************! 
      ! Calculo preditor elástico                                                !
      !**************************************************************************! 
      stresstrial = 0.0d0
      dstress     = MATMUL(dsdeEl(1:ncomp,1:ncomp),dstrain(1:ncomp))
      stresstrial = stressi + dstress
      !
      ! stresstrial = MATMUL(dsdeEl(1:ncomp,1:ncomp),
      !&              strain(1:ncomp)+dstrain(1:ncomp))
      !&              -epsEP(1:ncomp))
      !
      stress = stresstrial + sigi
      !
      !**************************************************************************! 
      ! Atualizando matriz constitutiva                                          !
      !**************************************************************************! 
      dsdePl = dsdeEl
      !
      !**************************************************************************! 
      ! Calcula função de escoamento                                             !
      !**************************************************************************!    
      CALL calcula_c(ci,cp,cr,eps1,eps2,eps3,epsEPeq,c)   ! coesão = f(epsEPeq)
      CALL invars(stress,ncomp,I1,J2,J3,p,q,theta,s)      ! Invariantes
      CALL yield(superficief,I1,J2,theta,c,fi,flag,f)     ! função de escoamento
      !
      !**************************************************************************! 
      ! Verifica o critério de escoamento                                        !
      !**************************************************************************! 
      IF(f.GT.0.0d0)THEN
          !
          ! Aplica o corretor plástico
          !
          !**************************************************************************! 
          ! Calcula dgds                                                             !
          !**************************************************************************!  
          dgds = 0.0d0
          CALL invars(stress,ncomp,I1,J2,J3,p,q,theta,s)      ! Invariantes do tensor de tensões 
          CALL g1g2g3(s,ncomp,J2,g1,g2,g3)                    ! Direção vetor de fluxo
          CALL c1c2c3(J2,theta,superficieg,psi,flag,c1,c2,c3) ! componentes vetor de fluxo
          dgds = c1*g1 + c2*g2 + c3*g3
          !
          !**************************************************************************! 
          ! Calcula h referente ao endurecimento e amolecimento                      !
          !**************************************************************************! 
          dhdq = 0.0d0
          CALL calcula_dfdc(superficief,flag,fi,dfdc)
          CALL calcula_dcde(ci,cp,cr,eps1,eps2,eps3,epsEPeq,dcde)
          dhdq = dfdc*dcde
          !
          !**************************************************************************! 
          ! Interações de NR local do modelo constitutivo                            !
          !**************************************************************************! 
          dlam = 0.0d0
          i = 0
          DO 
              !**************************************************************************! 
              ! Calculo de ddlamb                                                        !
              !**************************************************************************!
              ! Calcula dfds
              dfds = 0.0d0
              CALL invars(stress,ncomp,I1,J2,J3,p,q,theta,s)
              CALL g1g2g3(s,ncomp,J2,g1,g2,g3)
              CALL c1c2c3(J2,theta,superficief,fi,flag,c1,c2,c3)
              dfds    = c1*g1 + c2*g2 + c3*g3
              !
              ! Calcula dfdq
              dfdq = 0.0d0
              CALL calcula_dfdc(superficief,flag,fi,dfdc)
              CALL calcula_dcde(ci,cp,cr,eps1,eps2,eps3,epsEPeq,dcde)
              dfdq = dfdc
              !
              ! Calcula denominador de ddlamb
              ftDg    = DOT_PRODUCT(dfds,MATMUL(dsdeEl,dgds))
              !
              ! Calcula ddlamb
              ddlam   = f/(ftDg - dfdq*dhdq)
              !
              !**************************************************************************! 
              ! Calculo do corretor plástico                                             !
              !**************************************************************************!
              dstress = -ddlam*MATMUL(dsdeEl,dgds)
              !
              !**************************************************************************! 
              ! Incremento das tensões e do dlam                                         !
              !**************************************************************************!
              stress  = stress + dstress
              dlam    = dlam + ddlam
              !
              !**************************************************************************! 
              ! Calcula deformação plástica equivalente                                  !
              !**************************************************************************!
              depsEP = dlam*dgds
              CALL normatensor(depsEP,ncomp,depsEPeq)
              CALL calcula_Czao(flag,fi,Czao)
              epsEPeq = epsEPeq + Czao*depsEPeq
              !
              !**************************************************************************! 
              ! Verifica critério de escoamento                                          !
              !**************************************************************************!
              !CALL calcula_c(ci,cp,cr,eps1,eps2,eps3,epsEPeq,c)    ! coesão = f(epsEPeq)
              CALL invars(stress,ncomp,I1,J2,J3,p,q,theta,s)       ! Invariantes
              CALL yield(superficief,I1,J2,theta,c,fi,flag,f)      ! função de escoamento
              IF(f.LE.0.000001d0)EXIT
              !
              ! Caso atinja o número de iterações limites, faça a bisseção
              i = i + 1
              IF(i.EQ.nrmax)THEN
                  keycut = 1
                  RETURN
              ENDIF
          ENDDO
          !
          !**************************************************************************! 
          ! Atualizando o módulo constitutivo                                        !
          !**************************************************************************!       
          dfds_m(1,:) = dfds
          dgds_m(:,1) = dgds
          DgftD = MATMUL(MATMUL(MATMUL(dsdeEl,dgds_m),dfds_m),dsdeEl)
          dsdePl = dsdePl - DgftD/(ftDg - dfdq*dhdq)
          !
      ENDIF
      !
      !**************************************************************************! 
      ! Guardando deformações plásticas totais                                   !
      !**************************************************************************! 
      !
      
      !
      ! Calcula as deformações plásticas totais
      epsEP = epsEP + depsEP
      !
      ! Retorna as deformações inelasticas totais e equivalentes (só tem plasticas)
      epsPl = epsEP
      epseq = epsEPeq
      !
      ! Guarda deformação plástica equivalente
      ustatev(1) = epsEPeq
      !
      ! Guarda deformações plásticas
      CALL vmove(epsEP(1), ustatev(2), ncomp)
      !
      ! Guarda tensões
      CALL vmove(stress(1), ustatev(8), ncomp)
      !
      RETURN
      END
      !
c*deck,usermat3d    USERDISTRIB  parallel                                gal
c      subroutine usermat3d_EPVP (
c     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
c     &                   ldstep,isubst,keycut,
c     &                   nDirect,nShear,ncomp,nStatev,nProp,
c     &                   Time,dTime,Temp,dTemp,
c     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
c     &                   Strain,dStrain, epsPl, prop, coords, 
c     &                   var0, defGrad_t, defGrad,
c     &                   tsstif, epsZZ, cutFactor, 
c     &                   var1, var2, var3, var4, var5,
c     &                   var6, var7)
c#include "impcom.inc"  
c      !**************************************************************************!
c      !** SUBROTINA: USERMAT3D_EP                                              **!
c      !**                                                                      **!
c      !** Objetivo: atualiza as tensões, variáveis de estado e matriz          **!
c      !**           constitutiva para o modelo constitutivo elastoplastico     **!
c      !**                                                                      **!
c      !** Situação:                                                            **!
c      !**                                                                      **!
c      !**************************************************************************!         
c      ! 
c      !**************************************************************************! 
c      ! Declaração variáveis de entrada e saída da subrotina                     !
c      !**************************************************************************! 
c      INTEGER          
c     &                 matId, elemId,
c     &                 kDomIntPt, kLayer, kSectPt,
c     &                 ldstep,isubst,keycut,
c     &                 nDirect,nShear,ncomp,nStatev,nProp
c      DOUBLE PRECISION 
c     &                 Time,    dTime,   Temp,    dTemp,
c     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
c      DOUBLE PRECISION 
c     &                 stress  (ncomp  ), ustatev (nStatev),
c     &                 dsdePl  (ncomp,ncomp),
c     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
c     &                 epsPl   (ncomp  ), prop    (nProp  ), 
c     &                 coords  (3),
c     &                 defGrad (3,3),     defGrad_t(3,3),
c     &                 tsstif  (2)
c      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
c     &                 var6, var7
c      ! 
c      !**************************************************************************! 
c      ! Informações sobre as variáveis locais                                    !
c      !**************************************************************************! 
c      ! i,j           (int,sc,l)               contadores
c      ! mcomp         (int,sc,l)               total de componentes
c      ! stresstrial   (dp,ar(6),l)             tensão tentativa
c      ! dsdeEl        (dp,ar(6,6),l)           matriz constitutiva elastica
c      ! young         (dp,sc,l)                módulo de Young
c      ! posn          (dp,sc,l)                coeficiente de Poisson
c      ! 
c      !**************************************************************************! 
c      ! Declaração das variáveis locais                                          !
c      !**************************************************************************! 
c      INTEGER             i, j                ! contadores
c      INTEGER             mcomp               ! máximo de componentes
c      INTEGER             nrmax               ! máximo número de interações de NR
c      DOUBLE PRECISION    vPi                 ! Pi
c      PARAMETER           (mcomp = 6)         ! máximo de componentes
c      PARAMETER           (nrmax = 200)        ! máximo número de interações de NR
c      PARAMETER           (vPi = 3.14159265358979323846d0) ! máximo número de interações de NR
c      INTEGER             ncompgt             ! numero de componentes para o vetor de tensões iniciais
c      DOUBLE PRECISION    dstress(mcomp)      ! incremento de tensão
c      DOUBLE PRECISION    dsdeEl(mcomp,mcomp) ! matriz constitutiva elástica
c      DOUBLE PRECISION    dsdeElinv(mcomp,mcomp) ! inversa da matriz constitutiva elástica
c      DOUBLE PRECISION    sigi(mcomp)         ! tensões iniciais
c      !
c      DOUBLE PRECISION    young               ! modulo de Young
c      DOUBLE PRECISION    posn                ! coeficiente de Poisson
c      INTEGER             superficief         ! funçaõ de escoamento: 1- (MC ou TR), 2- (DP ou VM)
c      INTEGER             superficieg         ! função potencial: 1- (MC ou TR), 2- (DP ou VM)
c      INTEGER             flag                ! 1- DP coincide com as bordas externas MC 
c      DOUBLE PRECISION    c                   ! coesão
c      DOUBLE PRECISION    fi                  ! angulo de atrito
c      DOUBLE PRECISION    psi                 ! angulo de dilatância
c      DOUBLE PRECISION    ci,cp,cr            ! coesão inicial, pico e residual
c      DOUBLE PRECISION    eps1,eps2,eps3      ! deformações plásticas equivalentes limites
c      !
c      DOUBLE PRECISION    stresstrial(mcomp)  ! tensão tentativa
c      DOUBLE PRECISION    I1                  ! primeiro invariante das tensões
c      DOUBLE PRECISION    J2                  ! segundo invariante do desviador
c      DOUBLE PRECISION    J3                  ! terceiro invariante do desviador
c      DOUBLE PRECISION    p                   ! pressão hidrostática
c      DOUBLE PRECISION    q                   ! tensão equivalente
c      DOUBLE PRECISION    theta               ! angulo de Lode
c      DOUBLE PRECISION    s(mcomp)            ! tensor desviador
c      DOUBLE PRECISION    g1(mcomp),g2(mcomp),g3(mcomp)   ! direções do vetor de fluxo
c      DOUBLE PRECISION    c1,c2,c3            ! magnitudes do vetor de fluxo
c      DOUBLE PRECISION    dgds(mcomp)         ! vetor de fluxo
c      DOUBLE PRECISION    dfds(mcomp)         ! derivada de f em relação as tensões
c      DOUBLE PRECISION    f                   ! função de escoamento
c      DOUBLE PRECISION    ftDg                ! dfds*D*dgds
c      DOUBLE PRECISION    dlam                ! incremento da magnitude da deformação plástica
c      DOUBLE PRECISION    dlamk               ! incremento dlam pelo NR do corretor plástico
c      DOUBLE PRECISION    epsEP(mcomp)        ! deformação elastoplástica
c      DOUBLE PRECISION    depsEP(mcomp)       ! incremento de deformação plástica
c      DOUBLE PRECISION    depsEPeq            ! incremento de deformação plástica
c      DOUBLE PRECISION    epsEPeq             ! deformação plástica equivalente
c      DOUBLE PRECISION    czao                ! constante para o calculo da deformação plástica equivalente
c      DOUBLE PRECISION    DgftD(mcomp,mcomp)  ! Produto para calcular matriz algoritmica
c      DOUBLE PRECISION    dfds_m(1,mcomp)     ! Produto para calcular matriz algoritmica
c      DOUBLE PRECISION    dgds_m(mcomp,1)     ! Produto para calcular matriz algoritmica
c      DOUBLE PRECISION    dhdq,dfdq,dfdc,dcde
c      !
c      INTEGER             superficiefVP       ! funçaõ de escoamento: 1- (MC ou TR), 2- (DP ou VM)
c      INTEGER             superficiegVP       ! função potencial: 1- (MC ou TR), 2- (DP ou VM)
c      INTEGER             flagVP              ! 1- DP coincide com as bordas externas MC 
c      DOUBLE PRECISION    cVP                 ! coesão
c      DOUBLE PRECISION    fiVP                ! angulo de atrito
c      DOUBLE PRECISION    psiVP               ! angulo de dilatância
c      DOUBLE PRECISION    n,eta,f0            ! parâmetros do modelo de Perzyna
c      DOUBLE PRECISION    thetaVP             ! variável que controla o tipo de integração 0 - explicita, 1 - semimplicita
c      !
c      DOUBLE PRECISION    epsVP(mcomp)        ! deformação Viscoplástica
c      DOUBLE PRECISION    epsVPeq             ! deformação Viscoplástica equivalente
c      DOUBLE PRECISION    fvp                 ! superficie viscoplastica
c      DOUBLE PRECISION    PHI                 ! função de sobretensão
c      DOUBLE PRECISION    dgdsVP(mcomp)       ! vetor de fluxo viscoplástico
c      DOUBLE PRECISION    dPHIds(mcomp)       ! derivada de FI em relação as tensões
c      DOUBLE PRECISION    PHItDdeps           ! denominador do multiplicador viscoplástico
c      DOUBLE PRECISION    PHItDg              ! numerador do multiplicador viscoplástico
c      DOUBLE PRECISION    dlambdaVP          ! incremento do multiplicador viscoplástico
c      DOUBLE PRECISION    depsVP(mcomp)       ! incremento de deformação viscoplástica
c      DOUBLE PRECISION    dPHIds_m(1,mcomp)
c      DOUBLE PRECISION    dgdsVP_m(mcomp,1)
c      DOUBLE PRECISION    dgPHItD(mcomp,mcomp)
c      DOUBLE PRECISION    aux1
c      DOUBLE PRECISION    sigmap(mcomp)
c      !
c      DOUBLE PRECISION    depsEPVPeq
c      DOUBLE PRECISION    epsEPVPeq
c      DOUBLE PRECISION    depsPl(mcomp)
c      !
c      EXTERNAL            vzero, vmove, get_ElmData, get_ElmInfo
c      ! 
c      !**************************************************************************! 
c      ! Entrada de dados                                                         !
c      !**************************************************************************! 
c      keycut      = 0
c      young       = prop(2)
c      posn        = prop(3)
c      !
c      superficief = prop(4)
c      superficieg = prop(5)
c      flag        = prop(6)
c      fi          = prop(7)*vPi/180
c      psi         = prop(8)*vPi/180
c      ci          = prop(9)
c      cp          = prop(10)
c      cr          = prop(11)
c      eps1        = prop(12)
c      eps2        = prop(13)
c      eps3        = prop(14)
c      !
c      superficiefVP = prop(15)
c      superficiegVP = prop(16)
c      flagVP        = prop(17)
c      fiVP          = prop(18)*vPi/180
c      psiVP         = prop(19)*vPi/180
c      cVP           = prop(20)
c      n             = prop(21)
c      eta           = prop(22)
c      f0            = prop(23)
c      thetaVP       = prop(24)
c      ! 
c      !**************************************************************************! 
c      ! Coletando variáveis de estado                                            !
c      !**************************************************************************! 
c      epsEPeq = ustatev(1)
c      epsVPeq = ustatev(11)
c      CALL vmove(ustatev(2), epsEP(1), ncomp)
c      CALL vmove(ustatev(12), epsVP(1), ncomp)
c      ! 
c      !**************************************************************************! 
c      ! Calculo da matriz constitutiva                                           !
c      !**************************************************************************! 
c      dsdeEl   = 0.0d0
c      CALL MatrizD(young,posn,ncomp,dsdeEl)
c      ! 
c      !**************************************************************************! 
c      ! Calculo módulo de rigidez transversal para hourglass                     !
c      !**************************************************************************!  
c      tsstif(1) = 0.5d0*(young /(1.0d0+posn))
c      !
c      !**************************************************************************! 
c      ! Coletando tensões iniciais                                               !
c      !**************************************************************************! 
c      call get_ElmInfo('NCOMP', ncompgt)
c      call vzero(sigi(1),ncompgt)
c      call get_ElmData ('ISIG', elemId,kDomIntPt, ncompgt, sigi)  
c      !
c      !
c      !**************************************************************************! 
c      ! Calculo das deformações viscosas                                         !
c      !**************************************************************************!  
c      sigmap = 0.0d0
c      PHI = 0.0d0
c      CALL invars(stress,I1,J2,J3,p,q,theta,s)
c      CALL yield(superficiefVP,I1,J2,theta,cVP,fiVP,flagVP,fVP)
c      PHI = (fVP/f0)**n
c      IF(PHI.LE.0.0d0)PHI = 0.0d0
c      !
c      IF(PHI.GT.0.0d0)THEN
c          !**************************************************************************! 
c          ! Calculo dgdsvp                                                           !
c          !**************************************************************************! 
c          CALL g1g2g3(s,J2,g1,g2,g3)
c          CALL c1c2c3(J2,theta,superficiegVP,psiVP,flagVP,c1,c2,c3)
c          dgdsVP = c1*g1 + c2*g2 + c3*g3    
c          !
c          !**************************************************************************! 
c          ! Calculo dPHIds,PHItDg,PHItDdeps,aux1                                     !
c          !**************************************************************************! 
c          CALL c1c2c3(J2,theta,superficiefVP,psiVP,flagVP,c1,c2,c3)
c          dPHIds    = c1*g1 + c2*g2 + c3*g3
c          PHItDdeps = DOT_PRODUCT(dPHIds,MATMUL(dsdeEl,dstrain))
c          PHItDg    = DOT_PRODUCT(dPHIds,MATMUL(dsdeEl,dgdsVP))
c          aux1 = (eta/dtime + thetaVP*PHItDg)
c          !
c          !**************************************************************************! 
c          ! Calculo dlambda_vp                                                       !
c          !**************************************************************************!
c          dlambdaVP = (PHI + thetaVP*PHItDdeps)/aux1
c          depsVP = dlambdaVP*dgdsVP
c          !
c          ! Atualização do módulo constitutivo
c          dPHIds_m(1,:) = dPHIds
c          dgdsVP_m(:,1) = dgdsVP
c          DgPHItD = MATMUL(MATMUL
c     &              (MATMUL(dsdeEl,dgdsVP_m),dPHIds_m),dsdeEl)
c          dsdePl = dsdeEl - thetaVP*DgPHItD/aux1
c          !
c          ! Calculo da tensão inicial
c          sigmap = PHI*MATMUL(dsdeEl(1:ncomp,1:ncomp),
c     &            dgdsVP(1:ncomp))/aux1
c          !
c      ELSE
c          depsVP = 0.0d0  
c          dsdePl = dsdeEl
c      ENDIF
c      !
c      ! Guarda deformação viscoplástica equivalente
c      ustatev(11) = epsVPeq
c      !
c      ! Atualizando deformações viscplásticas
c      epsVP = epsVP + depsVP
c      !
c      ! Guarda deformações viscoplásticas
c      CALL vmove(epsVP(1), ustatev(12), ncomp)
c      !
c      !**************************************************************************! 
c      ! Calculo preditor elástico                                                !
c      !**************************************************************************! 
c      stresstrial = 0.0d0
c      stresstrial = MATMUL(dsdeEl(1:ncomp,1:ncomp),
c     &              strain(1:ncomp)+dstrain(1:ncomp)
c     &              -epsEP(1:ncomp)-epsVP(1:ncomp))
c      !
c      stress = stresstrial + sigi - sigmap
c      !
c      !**************************************************************************! 
c      ! Calcula dgds                                                             !
c      !**************************************************************************!  
c      ! Calcula os invariantes do tensor de tensões
c      CALL invars(stress,I1,J2,J3,p,q,theta,s)    
c      !
c      ! Calcula as direções do vetor de fluxo
c      CALL g1g2g3(s,J2,g1,g2,g3)
c      !
c      ! Calcula as magnitudes componentes do vetor de fluxo
c      CALL c1c2c3(J2,theta,superficieg,psi,flag,c1,c2,c3)
c      !
c      ! Calcula dgds
c      dgds = c1*g1 + c2*g2 + c3*g3
c      !
c      !**************************************************************************! 
c      ! Calcula função de escoamento                                             !
c      !**************************************************************************!  
c      CALL calcula_c(ci,cp,cr,eps1,eps2,eps3,epsEPeq,c)
c      CALL yield(superficief,I1,J2,theta,c,fi,flag,f)
c      !
c      !**************************************************************************! 
c      ! Verifica o critério de escoamento                                        !
c      !**************************************************************************! 
c      IF(f.GT.0.0d0)THEN
c          ! Aplica o corretor plástico
c          dlam = 0.0d0
c          i = 0
c          !
c          CALL calcula_dfdc(superficief,flag,fi,dfdc)
c          CALL calcula_dcde(ci,cp,cr,eps1,eps2,eps3,epsEPeq,dcde)
c          dhdq = dfdc*dcde
c          !
c          DO 
c              ! Calcula dfds
c              CALL invars(stress,I1,J2,J3,p,q,theta,s)
c              CALL g1g2g3(s,J2,g1,g2,g3)
c              CALL c1c2c3(J2,theta,superficief,fi,flag,c1,c2,c3)
c              dfds    = c1*g1 + c2*g2 + c3*g3
c              !
c              ! Calcula dlamda
c              ftDg    = DOT_PRODUCT(dfds,MATMUL(dsdeEl,dgds))
c              !
c              ! Calcula dfdq
c              CALL calcula_dfdc(superficief,flag,fi,dfdc)
c              CALL calcula_dcde(ci,cp,cr,eps1,eps2,eps3,epsEPeq,dcde)
c              dfdq = dfdc
c              !
c              dlamk   = f/(ftDg - dfdq*dhdq)
c              dstress = -dlamk*MATMUL(dsdeEl,dgds)
c              !
c              ! Atualizando variáveis internas
c              stress  = stress + dstress
c              dlam    = dlam + dlamk
c              !
c              ! calcula a deformação plástica equivalente
c              depsEP = dlam*dgds
c              CALL normatensor(depsEP,ncomp,depsEPeq)
c              CALL calcula_Czao(flag,fi,Czao)
c              epsEPeq = epsEPeq + Czao*depsEPeq
c              !
c              ! Testa o critério de escoamento
c              CALL invars(stress,I1,J2,J3,p,q,theta,s) 
c              CALL calcula_c(ci,cp,cr,eps1,eps2,eps3,epsEPeq,c)
c              CALL yield(superficief,I1,J2,theta,c,fi,flag,f)   
c              IF(f.LE.0.000001d0)EXIT
c              !
c              ! Caso atinja o número de iterações limites, faça a bisseção
c              i = i + 1
c              IF(i.EQ.nrmax)THEN
c                  keycut = 1
c                  RETURN
c              ENDIF
c          ENDDO
c          !
c          !**************************************************************************! 
c          ! Atualizando o módulo constitutivo                                        !
c          !**************************************************************************!       
c          dfds_m(1,:) = dfds
c          dgds_m(:,1) = dgds
c          DgftD = MATMUL(MATMUL(MATMUL(dsdeEl,dgds_m),dfds_m),dsdeEl)
c          dsdePl = dsdePl - DgftD/(ftDg - dfdq*dhdq)
c          !
c      ENDIF
c      !
c      ! calcula a deformação elastoplastica-viscoplastica equivalente
c      CALL normatensor(depsEP+depsVP,ncomp,depsEPVPeq)
c      CALL calcula_Czao(flag,fi,Czao)
c      epsEPVPeq = epsEPVPeq + Czao*depsEPVPpeq
c      !
c      ! Calcula as deformações elastoplasticas-viscoplasticas
c      depsPl = depsEP + depsVP
c      epsPl = epsPl + depsPl
c      !
c      ! Guarda deformação elastoplasticas-viscoplasticas equivalente
c      ustatev(1) = epsEPVPeq
c      !
c      ! Guarda deformações elastoplasticas-viscoplasticas
c      DO i=1,ncomp
c          ustatev(i+1) = epsPl(i)
c      ENDDO
c      !
c      RETURN
c      END
      !
      !
      !
      SUBROUTINE matrizD(E,Poisson,ncomp,D)
      !**************************************************************************!
      !** Função: matrizD                                                      **!
      !**                                                                      **!
      !** Objetivo: calcula a matriz consitutiva do material isotrópico        **!
      !**                                                                      **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      IMPLICIT NONE
      DOUBLE PRECISION E              ! módulo de elasticidade                  
      DOUBLE PRECISION Poisson        ! coeficiente de Poisson
      INTEGER ncomp                   ! numero de componentes
      DOUBLE PRECISION D(ncomp,ncomp) ! matriz constitutiva elástica isotrópica
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
	!
      IF(ncomp.GT.4)THEN
          D(ncomp-1,ncomp-1)=D(4,4)
	    D(ncomp,ncomp)=D(4,4)
      ENDIF
      !
      END SUBROUTINE MatrizD
      !
      SUBROUTINE invars(stress,ncomp,I1,J2,J3,p,q,theta,s)      
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
      DOUBLE PRECISION stress(ncomp)  ! tensões
      INTEGER          ncomp          ! numero de componentes
      DOUBLE PRECISION I1             ! primeiro invariante do tensor de tensões
      DOUBLE PRECISION J2,J3          ! primeiro e segundo invariante do desviador
      DOUBLE PRECISION s(6)           ! desviador
      DOUBLE PRECISION p,q            ! pressão hidrostática e tensão eq. de vm
      DOUBLE PRECISION theta          ! angulo de Lode
      DOUBLE PRECISION sine           ! variavel auxiliar
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
     &    stress(4)**2
      !
      IF(ncomp.GT.4)THEN
          J2 = J2 + stress(ncomp-1)**2+stress(ncomp)**2
      ENDIF
      !
      ! Calculo do p
      p = 1/3.0d0*I1
      !
      ! Calculo do desviador s
      s(1) = stress(1) - p
      s(2) = stress(2) - p
      s(3) = stress(3) - p
      s(4) = stress(4)
      IF(ncomp.GT.4)THEN
          s(ncomp-1) = stress(ncomp-1)
          s(ncomp) = stress(ncomp)
      ENDIF
      !
      ! Calculo do J3
      J3 = s(1)*s(2)*s(3)-s(3)*s(4)*s(4) 
      IF(ncomp.GT.4)THEN
          J3 = J3 -s(1)*s(ncomp-1)*s(ncomp-1)-s(2)*s(ncomp)*s(ncomp)+
     &        2.0d0*s(4)*s(ncomp-1)*s(ncomp)
      ENDIF
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
      SUBROUTINE g1g2g3(s,ncomp,J2,g1,g2,g3)
      !**************************************************************************!
      !** Subrotina: g1g2g3                                                    **!
      !**                                                                      **!
      !** Objetivo: calcula as direções do vetor de fluxo                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      IMPLICIT NONE
      DOUBLE PRECISION s(ncomp)                           ! desviador
      INTEGER          ncomp                              ! numero de componentes 
      DOUBLE PRECISION J2                                 ! segundo invariante do desviador
      DOUBLE PRECISION g1(ncomp), g2(ncomp), g3(ncomp)    ! direções do vetor de fluxo
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
      !
      ! Calculo do g2
      g2(1) = s(1)
      g2(2) = s(2)
      g2(3) = s(3)
      g2(4) = 2.0d0*s(4)
      IF(ncomp.GT.4)THEN
          g2(ncomp-1) = 2.0d0*s(ncomp-1)
          g2(ncomp) = 2.0d0*s(ncomp)
      ENDIF
      g2 = 1/(2*SQRT(J2))*g2
      !
      ! Calculo do g3
      g3(1) = s(2)*s(3) + J2/3.0d0
      g3(2) = s(1)*s(3) + J2/3.0d0
      g3(3) = s(1)*s(2) - s(4)**2 + J2/3.0d0
      g3(4) = 2.0d0*(-s(3)*s(4))
      IF(ncomp.GT.4)THEN
          g3(1) = g3(1) - s(ncomp-1)**2
          g3(2) = g3(2) - s(ncomp)**2 
          g3(4) = g3(4) + 2.0d0*s(ncomp-1)*s(ncomp)
          g3(ncomp-1) = 2.0d0*(s(ncomp)*s(4)-s(1)*s(ncomp-1))
          g3(ncomp) = 2.0d0*(s(4)*s(ncomp-1)-s(2)*s(ncomp)) 
      ENDIF      
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
      SUBROUTINE normatensor(tensor,ncomp,norma)      
      !**************************************************************************!
      !** Subrotina: normatensor                                               **!
      !**                                                                      **!
      !** Objetivo: calcula a norma de um tensor escrito em notação de Voigt   **!
      !**                                                                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!   
      IMPLICIT NONE
      INTEGER ncomp
      DOUBLE PRECISION tensor(ncomp)
      DOUBLE PRECISION norma
      !
      IF(ncomp.EQ.6)THEN
      norma = SQRT(2.0d0/6.0d0*
     &    ((tensor(1)-tensor(2))**2+(tensor(2)-tensor(3))**2+
     &    (tensor(3)-tensor(1))**2)+
     &    tensor(4)**2+tensor(5)**2+tensor(6)**2)
      ELSEIF(ncomp.EQ.4)THEN
          norma = SQRT(2.0d0/6.0d0*
     &    ((tensor(1)-tensor(2))**2+(tensor(2)-tensor(3))**2+
     &    (tensor(3)-tensor(1))**2)+
     &    tensor(4)**2)
      ELSE
      ENDIF
      !
      END SUBROUTINE
      !
      SUBROUTINE calcula_Czao(flag,fi,Czao)
      !**************************************************************************!
      !** Subrotina: Czao                                                      **!
      !**                                                                      **!
      !** Objetivo: calcula o C utilizado no calculo da deformação plástica    **!
      !**           efetiva                                                    **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************!           
      IMPLICIT NONE
      INTEGER             flag                ! 1- DP coincide com as bordas externas MC 
      DOUBLE PRECISION    fi                  ! angulo de atrito
      DOUBLE PRECISION    Czao                ! constante para o calculo da deformação plástica efetiva
      DOUBLE PRECISION    beta                ! constante referente a pressão hidrostática
      !
      ! Calcula o C para a deformação plástica equivalente
      IF(flag.EQ.1)THEN
          ! DP coincide com as bordas mais externas de MC
          beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0-SIN(fi)))
      ELSE
          ! DP coincide com as bordas mais internas de MC
          beta = 2.0d0*SIN(fi)/(SQRT(3.0d0)*(3.0d0+SIN(fi)))
      ENDIF          
      Czao = (beta+1.0d0/SQRT(3.0d0))/(SQRT(3*beta**2+1.0d0/2.0d0))
      END SUBROUTINE
      !
      SUBROUTINE calcula_dfdc(superficie,flag,fi,dfdc)
      !**************************************************************************!
      !** Subrotina: calcula_dfdc                                              **!
      !**                                                                      **!
      !** Objetivo: calcula dfdc                                               **!
      !**                                                                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************! 
      IMPLICIT NONE
      INTEGER             superficie          ! função potencial: 1- (MC ou TR), 2- (DP ou VM)
      INTEGER             flag                ! 1- DP coincide com as bordas externas MC 
      DOUBLE PRECISION    fi                  ! angulo de atrito
      DOUBLE PRECISION    dfdc                ! derivada de f em relação a c
      !
      SELECT CASE(superficie)
      CASE(1)
          dfdc = -COS(fi)
      CASE(2)
          SELECT CASE(flag)
          CASE(1)
              dfdc = -6*COS(fi)/(SQRT(3.0d0)*(3.0d0-SIN(fi)))
          CASE(2)
              dfdc = -6*COS(fi)/(SQRT(3.0d0)*(3.0d0+SIN(fi)))
          END SELECT
      END SELECT
      END SUBROUTINE
      !
      !
      SUBROUTINE calcula_dcde(ci,cp,cr,eps1,eps2,eps3,epsPleq,
     & dcde)
      !**************************************************************************!
      !** Subrotina: calcula_dfdepsPleq                                        **!
      !**                                                                      **!
      !** Objetivo: calcula df/depsPleq                                        **!
      !**                                                                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************! 
      IMPLICIT NONE
      DOUBLE PRECISION    ci                  ! coesão inicial
      DOUBLE PRECISION    cp                  ! coesão no pico
      DOUBLE PRECISION    cr                  ! coesão residual
      DOUBLE PRECISION    eps1                ! deformação plastica equivalente 1
      DOUBLE PRECISION    eps2                ! deformação plastica equivalente 2 
      DOUBLE PRECISION    eps3                ! deformação plastica equivalente 3
      DOUBLE PRECISION    epsPleq             ! deformação plástica equivalente
      DOUBLE PRECISION    dcde                ! dc/depsPleq
      !
      IF(epsPleq.LE.eps1)THEN
          dcde = (cp-ci)/(eps1)
      ELSEIF((epsPleq.GT.eps1).AND.(epsPleq.LE.eps2))THEN
          dcde = 0.0d0
      ELSEIF((epsPleq.GT.eps2).AND.(epsPleq.LE.eps3))THEN
          dcde = (cr-cp)/(eps3-eps2)
      ELSEIF(epsPleq.GT.eps3)THEN
          dcde = 0.0d0
      ENDIF
      !
      END SUBROUTINE
      !
      SUBROUTINE calcula_c(ci,cp,cr,eps1,eps2,eps3,epsPleq,
     & c)
      !**************************************************************************!
      !** Subrotina: calcul c                                                  **!
      !**                                                                      **!
      !** Objetivo: calcula c                                                  **!
      !**                                                                      **!
      !**                                                                      **!      
      !** Situação:                                                            **!
      !**                                                                      **!
      !**************************************************************************! 
      IMPLICIT NONE
      DOUBLE PRECISION    ci                  ! coesão inicial
      DOUBLE PRECISION    cp                  ! coesão no pico
      DOUBLE PRECISION    cr                  ! coesão residual
      DOUBLE PRECISION    eps1                ! deformação plastica equivalente 1
      DOUBLE PRECISION    eps2                ! deformação plastica equivalente 2 
      DOUBLE PRECISION    eps3                ! deformação plastica equivalente 3
      DOUBLE PRECISION    epsPleq             ! deformação plástica equivalente
      DOUBLE PRECISION    c                   ! coesão
      !
      IF(epsPleq.LE.eps1)THEN
          c = ci + (cp-ci)/(eps1)*epsPleq
      ELSEIF((epsPleq.GT.eps1).AND.(epsPleq.LE.eps2))THEN
          c = cp
      ELSEIF((epsPleq.GT.eps2).AND.(epsPleq.LE.eps3))THEN
          c = cp + (cr-cp)/(eps3-eps2)*(epsPleq - eps2)
      ELSEIF(epsPleq.GT.eps3)THEN
          c = cr
      ENDIF
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
