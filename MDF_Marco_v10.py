#!/usr/bin/env python
# *-* coding: utf-8 *-*
#https://www.wolframalpha.com/input/?i=u%3D%28%28%28x^4%29%2F12%29-%28%28x^3%29%2F6%29%29%2B%28%281%2F12%29%29x
#http://interactivepython.org/NWUbZ/LMdYZ/LpOMZ/courselib/static/pythonds/Graphs/graphintro.html

############################

import math
import numpy 
from scipy import linalg
import matplotlib.pyplot as plt



#Imprime uma lista em formato de matriz
def imprimeMatriz(matriz):
    print'\n'
    for i in range (len(matriz)):
        for j in range (len(matriz[0])):
            print matriz[i][j],
        print

#Imprime uma lista em formato de vetor
def imprimeVetor(vetor):
    print'\n'
    for i in range (len(vetor)):
        print vetor[i]

#Avalia se o problema pode ser resolvido
#Tc1 e Tc2 nao podem ser naturais simultaneamente
def verificacc(tc1,tc2):
    
    if tc1 != 1 and tc1 != 0:
        print 'Erro cc 1'
        
    if tc2 != 1 and tc2 != 0:
        print 'Erro cc 2'
        
    if tc1 == 1 and tc2 == 1:
        print 'Ambas as condicoes de contorno sao naturais => impossivel definir solucao'
        
    
#Avalia a entrada do numero de pontos
#Deve ser um valor inteiro >= 2
def verificanpts(npontos):
    if isinstance(npontos,int)== False:
        print 'Erro no numero de pontos (nao inteiro)'
    if npontos < 2:
        print 'Erro no numero de pontos (<2)'

#Montar a matriz de coeficientes da equação de laplace
#Via MDF
def MDFcoef1d(L,P,Propriedades,Trechos,npontos,tc1,vc1,tc2,vc2):

    h=(1.0*L)/(npontos-1) #h-dimensao dos trechos particionados
    
    VetorProp=[0]*npontos #cria um vetor com as propriedades dos meios, de tamanho
    for i in range (npontos): # npontos, a fim de facilitar a montagem do programa
        VetorProp[i] = Propriedades[P-1]
    b = 0
    for j in range (P-1):
        a = Trechos[j]/h    #obs: caso a discretizaçao coincida um ponto com a
        for k in range (b,int(a)+1): # transiçao de comportamento, a propriedade adotada
            VetorProp[k]=Propriedades[j] #sera a do trecho seguinte
        b = int(a)+1

    #imprimeVetor(VetorProp)

    matriz=[0]*npontos #linhas da matriz->cria uma lista com (npontos) elementos

    for i in range (npontos):
        matriz[i]=[0]*npontos #colunas da matriz->substitui cada elemento da lista por uma lista com (npontos) elementos

    if tc1==0: # 1a linha: condição de contorno essencial do ponto inicial       
        matriz[0][0]=1.0 
    elif tc1==1:# 1a linha: condição de contorno natural do ponto inicial
        matriz[0][0]=-1.0/h
        matriz[0][1]= 1.0/h


#preenche entre a linha [0] e a linha [(npontos)-1] da matriz de coeficientes
    if P == 1:
        for i in range(1,npontos-1):
            matriz[i][i+1]=1/h**2
            matriz[i][i-1]=1/h**2
            matriz[i][i]=(-2/h**2)

    if P!=1:
        for i in range (1,npontos-1):
            matriz[i][i+1]=(VetorProp[i]*(1/h**2))+((VetorProp[i+1]-VetorProp[i-1])/(4*h**2))
            matriz[i][i-1]=VetorProp[i]*(1/h**2)-((VetorProp[i+1]-VetorProp[i-1])/(4*h**2))
            matriz[i][i]=VetorProp[i]*(-2/h**2)


    if tc2==0: # última linha: condição de contorno essencial do ponto final       
        matriz[npontos-1][npontos-1]=1.0
    elif tc2==1:# ultima linha: condição de contorno natural do ponto final
        matriz[npontos-1][npontos-2]=-1.0/h
        matriz[npontos-1][npontos-1]= 1.0/h

    return matriz

#Monta o vetor C da equacao matricial (MDFcoef1d)X = C
#Sendo X o vetor com os valores de x onde se deseja definir u(x):

#Problema de Laplace
def vetorlaplace(npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos

    vetor[0]=vc1
    vetor[npontos-1]=vc2

    return vetor

#Problema de Poisson para f(x) = x² - x
def vetorpoisson(k2,L,npontos,tc1,vc1,tc2,vc2):
    vetor=[0]*npontos
    h=(1.0*L)/(npontos-1)
    vetor[0]=vc1

    if k2 == 1:
        for i in range (1,npontos-1):
            x=i*h
            funcao=(x**2)-x
            vetor[i]=funcao
    if k2 == 2:
        for i in range (1,npontos-1):
            x=i*h
            funcao=x
            vetor[i]=funcao
    if k2 == 3:
        for i in range (1,npontos-1):
            x=i*h
            funcao=x*(x-1)*(x-0.5)
            vetor[i]=funcao
    if k2 == 4:
        for i in range (1,npontos-1):
            x=i*h
            funcao=math.exp(x)
            vetor[i]=funcao
    if k2 == 5:
        for i in range (1,npontos-1):
            x=i*h
            funcao=math.sin(x) + math.cos(x)
            vetor[i]=funcao
    if k2 == 6:
        for i in range (1,npontos-1):
            x=i*h
            funcao=x
            vetor[i]=funcao
    if k2 == 7:
        for i in range (1,npontos-1):
            x=i*h
            if x < 0.5*L:
                funcao=0
            if x >= 0.5*L:
                funcao=1
            vetor[i]=funcao
        
    vetor[npontos-1]=vc2

    return vetor

#Montar o vetor resposta U(x) para os pontos x desejados
#Laplace solucionado via analitica
def exatalaplace1d(L,npontos,tc1,vc1,tc2,vc2):  
    h=1.0*L/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
         for i in range (npontos):
            x=i*h
            exata[i]=(vc2-vc1)*x+vc1

    if tc1==0 and tc2==1:
         for i in range (npontos):
            x=i*h
            exata[i]=vc2*x+vc1

    if tc1==1 and tc2==0:
         for i in range (npontos):
            x=i*h
            exata[i]= vc1*x+(vc2-vc1)

    return exata

#Calculo da solucao analitica do problema de Poisson#########################
#Montar o vetor resposta U(x) para os pontos x desejados
def exatapoisson1d(k2,Propriedades,npontos,L,tc1,vc1,tc2,vc2):  
    h=1.0/(npontos-1)
    exata=[0]*npontos
       
    if tc1==0 and tc2==0:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+((1.0/12.0)+vc2-vc1)*x + vc1
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)-((1.0/6.0) + vc1 - vc2)*x + vc1
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) -(vc1-vc2+(1.0/120.0))*x + vc1
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (math.e - 1 + vc1 - vc2)*x - 1+vc1
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (1+vc1-vc2-math.sin(1)-math.cos(1))*x +1+vc1
                exata[i]=u
        if k2==6:
            h=1.0*L/(npontos-1)
            exata=[0]*npontos
            K1=Propriedades[0]
            K2=Propriedades[1]
            for i in range (npontos):
                x=i*h
                c3 = 1.0*vc1
                c1 = ((((L**2)/24.0)*((1.0/K2)-(1.0/K1)))+((2.0/L)*((vc2-vc1))-((L**2)/(3.0*K2))))/(1.0+(K1/(1.0*K2)))
                c2 = (1.0)*c1*K1/K2
                c4 = ((1.0)*vc2)-((L**3)/(6.0*K2))-(L*c2)
                if x<0.5*L:
                    u=((x**3)/(6.0*K1))+(c1*x)+c3
                if x>=0.5*L:
                    u=((x**3)/(6.0*K2))+(c2*x)+c4
                exata[i]=u
        if k2==7:
            h=1.0*L/(npontos-1)
            exata=[0]*npontos
            for i in range (npontos):
                x=i*h
                c3 = ((-5.0*L)/8.0)+((vc2-vc1)/(L*1.0))
                c1 = (L/2.0)+c3
                c4 = (vc2*1.0)-((L**2)/2.0)-(c3*L)
                c2 = vc1*1.0
                if x<0.5*L:
                    u=c1*x+c2
                if x>=0.5*L:
                    u=((x**2)/2.0)+(c3*x)+c4
                exata[i]=u
                    
            
    if tc1==0 and tc2==1:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+(vc2 + (1/6.0))*x + vc1
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)-(0.5-vc2)*x + vc1
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) -((1.0/120.0)-vc2)*x + vc1
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (math.e-vc2)*x - 1+vc1
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (math.sin(1)-math.cos(1)-vc2)*x +1+vc1
                exata[i]=u

    if tc1==1 and tc2==0:
        if k2==1:
            for i in range (npontos):
                x=i*h
                u=(((x**4)/12.0)-((x**3)/6.0))+(vc1)*x + ((1/12.0)+vc2-vc1)
                exata[i]=u
        if k2==2:
            for i in range (npontos):
                x=i*h
                u=((x**3)/6.0)+(vc1)*x -((1.0/6.0) + vc1 - vc2)
                exata[i]=u
        if k2==3:
            for i in range (npontos):
                x=i*h
                u=(((x**5)/20)-((x**4)/8)+((x**3)/12)) + (vc1)*x - (vc1-vc2+(1.0/120.0))
                exata[i]=u
        if k2==4:
            for i in range (npontos):
                x=i*h
                u=(math.exp(x)) - (1-vc1)*x - (math.e - 1 + vc1 - vc2)
                exata[i]=u
        if k2==5:
            for i in range (npontos):
                x=i*h
                u=-(math.sin(x)+math.cos(x)) - (1+vc1)*x - (1+vc1-vc2-math.sin(1)-math.cos(1))
                exata[i]=u

    return exata

def erroabsoluto(solnumerica,solexata,npontos):
    erro=[0]*npontos
    for i in range (npontos):
        erro[i]=solnumerica[i]-solexata[i]
    return erro
#Calculo da solucao analitica do problema de Poisson#########################

#Avaliacao de erro###########################################################
def calculaarea(solucao,npontos):
    area = 0
    h = (1.0*L)/(npontos-1)
    for i in range (npontos-1):
        area = ((solucao[i+1]+solucao[i])*.5*h) + area
    return area

def calculaerro(solnumerica,solnumanterior,npontos,refinamento,erroatual):
    if refinamento == 0:
        erroanterior = 999999
    if refinamento != 0:
        erroanterior = erroatual
    areaatual = calculaarea(solnumerica,npontos)
    erroatual = areaatual - areaanterior
    if erroatual < 0:
        erroatual = (-1.0)*erroatual
    solnumanterior = solnumerica
    avaliacao = erroatual/areaanterior
    if avaliacao < 0:
        avaliacao = (-1.0)*avaliacao
    return erroanterior, erroatual, solnumanterior, avaliacao


def avaliaerro(erroanterior,erroatual,npontos,refinamento,parada,avaliacao):
    diferro = erroatual - erroanterior
    if diferro<0:
        diferro = (-1.0)*diferro
    if avaliacao <= erroadm:
        parada = 1
    refinamento = refinamento + 1
    return npontos, refinamento, parada, diferro
#Avaliacao de erro###########################################################    
    

##############################################################################
#Estrutura principal do programa:
##############################################################################
k1=input('Digite "1" para resolver Laplace; "2" para Poisson:')

if k1==2:
    k2=input('Escolha uma das funçoes da lista:\n'
             'a)Para o caso homogeneo:\n'#caso compare o erro pela analitica, faça L=1
             '(1) => f(x) = x(x-1)\n'
             '(2) => f(x) = x\n'
             '(3) => f(x) = x(x-1)(x-0,5)\n'
             '(4) => f(x) = e**x\n'
             '(5) => f(x) = senx + cosx\n'
             'b)Para o caso heterogeneo:\n'
             '(6)k1!=k2, f(x)=x, u1 e u2 prescritos, transicao em L/2\n'
             '(7)k1=k2, f1(x)=0, f2(x)=1, u1 e u2 prescritos, transicao em L/2\n')

L=input('Tamanho do dominio do problema:')


P=input('O meio é composto por quantos materiais?')
Propriedades = [0]*P
for i in range (P):
    print 'Propriedade',i+1,':',
    Propriedades[i]= input('')
Trechos = [0]*(P-1)
for i in range (P-1):
    print 'Ponto de transição',i+1,':',
    Trechos[i]=input('')
    
npontos=input('Pontos para iniciar a discretizacao:')

verificanpts(npontos)

tc1=input('Tipo de condicao de contorno da esquerda:')
vc1=input('Valor de condicao de contorno da esquerda:')

tc2=input('Tipo de condicao de contorno da direita:')
vc2=input('Valor de condicao de contorno da direita:')

verificacc(tc1,tc2)
#################################################
#Iteracao zero
#################################################
matrizcoef=MDFcoef1d(L,P,Propriedades,Trechos,npontos,tc1,vc1,tc2,vc2)

if k1==1:
    vetorMDF=vetorlaplace(npontos,tc1,vc1,tc2,vc2)
    solexata=[0]*npontos#apenas para funcionamento do programa
elif k1==2:
    vetorMDF=vetorpoisson(k2,L,npontos,tc1,vc1,tc2,vc2)
    solexata=exatapoisson1d(k2,Propriedades,npontos,L,tc1,vc1,tc2,vc2)

#imprimeMatriz(matrizcoef)

solnumanterior=linalg.solve(matrizcoef,vetorMDF)
erroreal=erroabsoluto(solnumanterior,solexata,npontos)
areaanterior = calculaarea(solnumanterior,npontos)

#################################################
#Escolha sobre comparaçao analitica e numerica ou pela sucessao de refinamentos
#################################################
print '\n'
k3=input('Deseja avaliar o erro com base:\n'
         '(1)na solucao analitica;\n'
         '(2)pelo metodos de refinamentos.\n')

if k3==1:
    abscissas=numpy.linspace(0,L,npontos)
    imprimeVetor(solnumanterior)
    print 'Soluçao Numerica'
    imprimeVetor(solexata)
    print 'Soluçao Exata\n'
    imprimeVetor(erroreal)
    print 'Erro absoluto\n'
    plt.plot(abscissas,solnumanterior,'bo-', label='Numerica')   #gráfico soluçao numerica
    plt.plot(abscissas,solexata,'xr-', label='Exata')   #gráfico solucao exata
    plt.legend(loc='upper left')
    #Identificaçao dos eixos
    axes = plt.gca()
    axes.set_xlabel('x')
    axes.set_ylabel('u(x)')
    plt.show()    #apresenta os gráficos plotados
    
if k3==2:

    refinamento = 0

    parada = 0

    passo = 10

    erroadm = .5

    erroatual = 0#valor atribuido apenas para definir erroatual como variavel global

    while parada == 0:

        npontos = npontos + passo

        matrizcoef=MDFcoef1d(L,P,Propriedades,Trechos,npontos,tc1,vc1,tc2,vc2)

        #imprimeMatriz(matrizcoef)
        #print 'Matriz de Coeficientes MDF'

        if k1==1:
            vetorMDF=vetorlaplace(npontos,tc1,vc1,tc2,vc2)

        elif k1==2:
            vetorMDF=vetorpoisson(k2,L,npontos,tc1,vc1,tc2,vc2)

        #imprimeVetor(vetorMDF)
        #print 'Vetor MDF'

        solnumerica=linalg.solve(matrizcoef,vetorMDF)

        #imprimeVetor(solnumerica)
        #print 'Soluçao Numerica'

        erroanterior, erroatual, solnumanterior, avaliacao = calculaerro(solnumerica,solnumanterior,npontos,refinamento,erroatual)

        npontos, refinamento, parada, diferro = avaliaerro(erroanterior,erroatual,npontos,refinamento,parada,avaliacao)
    
    print '\n'
    print  'Erro atual:', erroatual
    print '\n'
    print 'Erro anterior:', erroanterior
    print '\n'
    print 'Quantidade de pontos para discretizacao final:', npontos
    print '\n'
    print 'Quantidade de refinamentos:', refinamento
    print '\n'
    print 'Diferenca entre valores de erro:', diferro
    print '\n'

    abscissas=numpy.linspace(0,L,npontos)

    #imprimeVetor(abscissas)
    #print 'Pontos discretizados'

    plt.plot(abscissas,solnumerica,'bo-', label= 'Numerica')

    plt.legend(loc='upper left')

    #Identificaçao dos eixos
    axes = plt.gca()
    axes.set_xlabel('x')
    axes.set_ylabel('u(x)')


    plt.show()    #apresenta os gráficos plotados
