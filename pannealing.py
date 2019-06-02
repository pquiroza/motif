import random
import copy
import decimal
import sys
import uuid
import hashlib


from multiprocessing import Process, Queue
from collections import Counter

datos = []
datos2 = []
datos3 = []
candidatos = []
cadn = []
limite = decimal.Decimal("0.95")

import math

class Secuencia:
    def __init__(self,nombre,seq):
        self.nombre = nombre
        self.seq = seq


class Kmer:
    def __init__(self,idgenome,adn,posicioninicial,largo):
        self.idgenome=idgenome
        self.adn=adn
        self.posicioninicial = posicioninicial
        self.largo = largo



class matrix:
    def __init__(self, kmers,fitness,fitness2,motif,code):
        self.kmers = kmers
        self.fitness = fitness
        self.fitness2 = fitness2
        self.motif = motif
        self.code = code

    def setCode(self,code):
        self.code = code

    def setKmer(self,posicion,kmer):
        self.kmers[posicion] = kmer

    def setFitness(self,fitness):
        self.fitness = fitness

    def setFitness2(self,fitness2):
        self.fitness2 = fitness2

    def setMotif(self,motif):
        self.motif = motif



class Candidato:
    def __init__(self,matrix,adn):
        self.matrix = matrix
        self.adn = adn

def HammingDistance(kmer1,kmer2):
    contador = 0
    for i in range(len(kmer1)):
        if (kmer1[i]!=kmer2[i]):
            contador += 1
    return contador


def HammingMatrix(matriz,datos):
    contador = 0
    totales = []
    for i in matriz.kmers:
        hamming = []

        for j in matriz.kmers:
            distancia =  HammingDistance(i.adn,j.adn)
            if (distancia==0 and i!=j):
                contador = contador + 1

            hamming.append(distancia)
        suma = sum(hamming)

        totales.append(suma)


    if (contador/2 > 2):
        adn = valorrepetido(matriz)
        valor = 0

        for cd in cadn:
            if (cd == adn):
                valor = 1


        if (valor==0):
            candidatos.append(matriz)
            cadn.append(adn)


        #escribemejor("candidatos2.fts",matriz,len(matriz.kmers[0].adn),"CANDIDATO")

    total = sum(totales)

    return total

def HammingMatrixr(matriz,datos):
    contador = 0
    totales = []


    for i in range(len(matriz.kmers)):
        hamming = []

        for j in range(i,len(matriz.kmers)):
            distancia =  HammingDistance(matriz.kmers[i].adn,matriz.kmers[j].adn)
            if (distancia==0 and i!=j):
                contador = contador + 1

            hamming.append(distancia)
        suma = sum(hamming)

        totales.append(suma)


    if (contador/2 > 2):
        adn = valorrepetido(matriz)
        valor = 0

        for cd in cadn:
            if (cd == adn):
                valor = 1


        if (valor==0):
            candidatos.append(matriz)
            cadn.append(adn)


        #escribemejor("candidatos2.fts",matriz,len(matriz.kmers[0].adn),"CANDIDATO")

    total = sum(totales)

    return total

def HammingMatrix2(matriz):

    totales = []
    for i in matriz.kmers:
        hamming = []
        for j in matriz.kmers:


            distancia =  HammingDistance(i.adn,j.adn)
            hamming.append(distancia)
        suma = sum(hamming)

        totales.append(suma)


    return totales

def getPalabra(datos,start,posicioninicial,largo):
    genome = datos[start]
    palabra = ""
    for i in range(posicioninicial,posicioninicial+largo):
        palabra = palabra+genome[i]

    return palabra

def generaVecino(matrix,datos,largo):
    vecino = copy.deepcopy(matrix)
    cambios = random.randint(1,len(vecino.kmers)-1)


    start = random.randint(0,len(vecino.kmers)-1)



    kmer = Kmer(start,vecino.kmers[start].adn,vecino.kmers[start].posicioninicial,vecino.kmers[start].largo)
    posicion = random.randint(0,len(kmer.adn))
    kmer.posicioninicial = posicion
    kmer.adn = getPalabra(datos,start,kmer.posicioninicial,kmer.largo)

    vecino.setKmer(start,kmer)


    return vecino

def escribeFitness(fitness,lkmer,archivo):
    f = open(archivo,'a')
    f.write(str(fitness))
    f.write('\n')
    f.close

def escribefinal(archivo,candidato,lkmer,mensaje):
    f=open(archivo,'a')
    f.write("RESULTADO "+str(lkmer)+" "+mensaje)
    f.write('\n')
    for i in candidato.kmers:
        f.write(i)
        f.write('\n')

    f.write("---------")
    f.write('\n')
    f.close

def escribemejor(archivo,mejor,lkmer,mensaje):
    f=open(archivo,'a')
    f.write("Motifs "+str(lkmer)+" "+mensaje)
    f.write('\n')
    for i in mejor.kmers:

        f.write(i.adn+" PI " +str(i.posicioninicial))
        f.write('\n')

    f.write("HM " + str(mejor.fitness))
    f.write('\n')
    f.write("MF " + str(mejor.fitness2))
    f.write('\n')
    f.write("MOTIF " + mejor.motif)
    f.write('\n')
    f.write("---------")
    f.write('\n')
    f.close



def generaVecino2(datos,largo):
    kmers = []
    matriz = matrix(kmers,0,0,"")
    identificador = 0
    for i in datos:
        inicio = random.randint(0,len(i)-largo)
        palabra = ""
        for j in range(inicio,inicio+largo):

            palabra = palabra + i[j]


        kmer = Kmer(identificador,palabra,inicio,largo)
        identificador += 1
        matriz.kmers.append(kmer)

    #for i in matriz.kmers:

    return matriz

def generaMatrizInicial(datos,largo):
    kmers = []
    matriz = matrix(kmers,0,0,"","")

    identificador = 0
    for i in datos:
        inicio = random.randint(0,len(i)-largo)
        palabra = ""
        for j in range(inicio,inicio+largo):

            palabra = palabra + i[j]

        kmer = Kmer(identificador,palabra,inicio,largo)
        identificador += 1
        matriz.kmers.append(kmer)



    return matriz


def fitnessAnnealing(matriz):
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    fitness = 0
    palabra = ""
    nkmers = len(matriz.kmers)
    lkmers = len(matriz.kmers[0].adn)
    for k in range(len(matriz.kmers[0].adn)):
        pal = ""
        counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}
        for i in matriz.kmers:
            counts[i.adn[k]]+=1
            pal = pal + i.adn[k]
        valor = counts.values()


        maximo = max(valor)

        if (maximo == nkmers):
            maximo = maximo + lkmers
        fitness = fitness + maximo

        ceros = list(counts.values()).count(0)

        fitness = fitness + ceros



    propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal((len(datos)+lkmers+20)*lkmers)))
    matriz.setFitness(propor)

    return propor



def parallelFitness(mmatriz,q):
    localr = []
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    for matriz in mmatriz:
        fitness = 0
        palabra = ""
        nkmers = len(matriz.kmers)
        lkmers = len(matriz.kmers[0].adn)
        for k in range(len(matriz.kmers[0].adn)):
            pal = ""
            counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}
            for i in matriz.kmers:
                counts[i.adn[k]]+=1
                pal = pal + i.adn

            ho = hashlib.md5(pal.encode())
            valor = counts.values()


            maximo = max(valor)

            if (maximo == nkmers):
                maximo = maximo + lkmers
            fitness = fitness + maximo

            ceros = list(counts.values()).count(0)

            fitness = fitness + ceros



        propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal((len(datos)+lkmers+20)*lkmers)))
        matriz.setFitness(propor)
        matriz.setCode(ho.hexdigest())


        localr.append(matriz)

    q.put(localr)



def getFitness(matriz):
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]

    fitness = 0
    palabra = ""
    nkmers = len(matriz.kmers)
    for k in range (len(matriz.kmers[0].adn)):

        counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}
        for i in matriz.kmers:

            counts[i.adn[k]]+=1
        valor =  counts.values()
        maximo =  max(valor)


        #posicion = valor.index(maximo)


        fitness = fitness + maximo

        #palabra = palabra + alfabeto[posicion]





    return fitness, palabra

def showMatriz(matriz):
    for i in matriz.kmers:
        print (i.adn)
    print (matriz.fitness)


def pwm(matrix,lkmer):
    score = 0
    pwm = []
    final = []
    probfinal = []
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y"]
    nkmers = len(matrix.kmers)
    for k in range (len(matrix.kmers[0].adn)):
        counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
        for i in matrix.kmers:
            counts[i.adn[k]]+=1
        #print counts
        #print counts.values()
        pwm.append(counts.values())
        #print "-" * 10
    valor =  counts.values()
    #for v in matrix.kmers:
        #   print v.adn
    #print pwm
    for i in range (len(pwm[0])):
        parcial = []
        for j in pwm:
            parcial.append(j[i])
        final.append(parcial)

    for fila in final:
        prob = []
        for j in fila:
            valor = decimal.Decimal(j)/decimal.Decimal(lkmer)
            razon = decimal.Decimal(1)/decimal.Decimal(len(alfabeto))


            if (valor!=0):

                total = math.log(decimal.Decimal(valor)/decimal.Decimal(razon))
            else:
                total = 0
            prob.append(str(total))
        probfinal.append(prob)
    for p in probfinal:
        print (p)
        print ("\n")
    return probfinal


def valorpwm(candi,pwm):
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y"]
    counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}
    palabra = candi.motif

    score = 0.0
    for i in range (len(candi.motif)):
        print (candi.motif[i])
        indice =   alfabeto.index(candi.motif[i])
        print (indice)

        score = pwm[indice][i]
        print (score)


def valorrepetido(matrix):
    valores = []
    for i in matrix.kmers:
        cont = 0
        for j in matrix.kmers:
            if i.adn==j.adn:
                cont = cont +1
        valores.append(cont)
        maximo = max(valores)
        valor = valores.index(maximo)
    #print "CANDIDATO " +  matrix.kmers[valor].adn
    return matrix.kmers[valor].adn

def procesaCandidatos(candidatos,datos,lkmer,umbral,salida):




    for c in candidatos:

        candi = copy.deepcopy(c)

        adn = valorrepetido(candi)




        for i in range (len(datos)):
            distanciamenor = lkmer+1
            palabrafinal = ""
            pos = 0
            for d in range(len(datos[i])-lkmer+1):
                palabra = ""


                for j in range (d,d+lkmer):
                    palabra = palabra + datos[i][j]

                distancia = HammingDistance(adn,palabra)


                if (distancia<distanciamenor):

                    pos = d
                    palabrafinal = palabra
                    distanciamenor = distancia

            kmerfinal = Kmer(i,palabrafinal,pos,lkmer)
            candi.setKmer(i,kmerfinal)




        fitness = HammingMatrixr(candi,datos)

        ft, pl = getFitness(candi)



        propor = decimal.Decimal(decimal.Decimal(ft)/(decimal.Decimal(lkmer*len(datos))))

        candi.setFitness2(propor)
        candi.setMotif(pl)




        candi.setFitness(fitness)

        limite = decimal.Decimal(umbral)
        if (propor>limite):
            escribemejor(salida,candi,lkmer,"Final ")





#>sp|P15783|PSPC_BOVIN Pulmonary surfactant-associated protein C OS=Bos taurus GN=SFTPC PE=1 SV=3
def cargaSecuencias(archivo):
    datos = []
    f = open(archivo,'r')
    secuencia = ""
    for line in f:
        p = line
        if p[0]==">":
                if (secuencia!=""):
                    datos.append(secuencia)
                    secuencia = ""
                else:
                    sp = p.split("|",2)
                    ns = sp[2].split(" ",2)
                    nombre = ns[0]

                    secuencia = ""

        else:
            l = line.split("\n",1)

            secuencia = secuencia + l[0]
    datos.append(secuencia)



    return datos

def ejecutar(ciclo,tinicial,lkmer,datos,salida,umbral):
    fitnessmejor = 10000
    #limite = (lkmer * len(datos)) - (lkmer*1.5)

    m=generaMatrizInicial(datos,lkmer)
    mejor = generaMatrizInicial(datos,lkmer)
    mejor.setFitness(1000)
    #fitnessinicial,palabrainicial = getFitness(m)
    fitnessinicial = HammingMatrix(m,datos)

    while (tinicial>0.2):

        #fitnessinicial,palabrainicial = getFitness(m)
        fitnessinicial = HammingMatrix(m,datos)

        #print "Fitness : " + str(fitnessinicial) + " " + palabrainicial  + " " + str(tinicial) + " fitness mejor " + str(fitnessmejor)
        print (str(tinicial) +" "+str(lkmer))
        for i in range (ciclo):


            vecino = generaVecino2(datos,lkmer)





            fitnessnuevo = HammingMatrix(vecino,datos)
            vecino.setFitness(fitnessnuevo)

            fitnessactual = HammingMatrix(m,datos)
            m.setFitness(fitnessactual)
            escribeFitness(fitnessactual,lkmer,"results/annealing")


            if (m.fitness<mejor.fitness):
                print ("MEJOR a " + str(m.fitness))
                mejor = copy.deepcopy(m)
                if (mejor.fitness<100):
                    escribemejor("Mejor.fts",mejor,lkmer,"Mejor")

            diferencia = vecino.fitness-m.fitness

            if (diferencia<0):

                m = copy.deepcopy(vecino)
            else:
                rnd = random.uniform(0.0,1.0)
                valor = (-diferencia) / decimal.Decimal(str(tinicial))
                exp = math.exp(valor)

                if  (rnd<exp):

                    m=copy.deepcopy(vecino)

            #print "-"*10

        tinicial = tinicial*0.95;
    print ("Numero Candidatos " +  str(len(candidatos)))





    procesaCandidatos(candidatos,datos,lkmer,umbral,salida)














def annealing2(ciclo,tinicial,lkmer,datos,salida):
    m = generaMatrizInicial(datos,lkmer)
    mejor = generaMatrizInicial(datos,lkmer)
    mejor.setFitness(1000)
    #fitnessinicial,palabrainicial = getFitness(m)


    while (tinicial>0.1):
        fitnessinicial = fitnessAnnealing(m)
        m.setFitness(fitnessinicial)
        print (str(tinicial) +" "+str(lkmer))
        for i in range(ciclo):
            fitnessactual = fitnessAnnealing(m)
            m.setFitness(fitnessactual)
            #showMatriz(m)
            vecino = copy.deepcopy(m)
            vecino = mutacion(vecino,datos)

            fitnessnuevo = fitnessAnnealing(vecino)
            vecino.setFitness(fitnessnuevo)
            #showMatriz(vecino)


            escribeFitness(fitnessactual,lkmer,"results/annealing")

            diferencia = decimal.Decimal(m.fitness)-decimal.Decimal(vecino.fitness)
            #print fitnessactual,fitnessnuevo,diferencia
            #print(m.fitness,vecino.fitness)
            #print(diferencia)


            #print diferencia,tinicial
            if (diferencia<0):
                #print("eligemejor")
                m = copy.deepcopy(vecino)
            else:
                rnd = random.uniform(0.0,1.0)
                valor = decimal.Decimal(-diferencia) / decimal.Decimal(str(tinicial))
                exp = math.exp(valor)
                print(tinicial,rnd,exp)
                if  (rnd<exp):
                    #print("eligepeor")
                    m=copy.deepcopy(vecino)


        tinicial = tinicial * 0.9
    #print ("Canditados " + str(len(candidatos)))
    #procesaCandidatos(candidatos,datos,lkmer,umbral,salida)




def cruzamiento(m1,m2):

    hijo1 = copy.deepcopy(m1)
    hijo2 = copy.deepcopy(m2)

    cruce1 = random.randint(1,len(m1.kmers)-1)



    for i in range(cruce1,len(m1.kmers)):
        kmerp = m1.kmers[i]
        hijo1.setKmer(i,m2.kmers[i])
        hijo2.setKmer(i,kmerp)


    return hijo1,hijo2



def escribeindividuo(individuo,generacion,archivo):
    f=open(archivo,'a')
    f.write("Motifs "+str(generacion))
    f.write('\n')
    for i in individuo.kmers:

        f.write(i.adn +" " +str(i.posicioninicial))
        f.write('\n')

    f.write(str(individuo.fitness))
    f.write('\n')
    f.close

def newseleccion(population,umbral,media):
    mejores = []
    print(media)
    valor = decimal.Decimal(media)*decimal.Decimal(1.1)
    for i in population:
        if i.fitness>=(valor):

            mejores.append(i)


    return mejores


def parallelSeleccion(poplocal,q):
    ganadores = []
    distintos = {}
    largo = len(poplocal)
    while(len(poplocal)>2):
        participante1 = random.randint(0,len(poplocal)-2)
        participante2 = random.randint(0,len(poplocal)-2)

        if (poplocal[participante1].fitness >= poplocal[participante2].fitness):
            if (poplocal[participante1].code not in distintos):
                ganadores.append(poplocal[participante1])
                distintos[poplocal[participante1].code]=1



        else:
            if (poplocal[participante2].code  not in distintos):
                ganadores.append(poplocal[participante2])
                distintos[poplocal[participante1].code]=1

        poplocal.pop(participante1)
        poplocal.pop(participante2)

    q.put(ganadores)

def seleccion(population,umbral):

    mejores = []
    for j in range(umbral):
        indice = 0
        fm = 0
        valor = 0

        for i in population:

            if i.fitness>fm:
                mejor = i
                fm = i.fitness
                valor = indice

            indice = indice + 1
        #print str(fm) + " " + str(j)
        mejores.append(mejor)
        if (len(population)>0):
            population.pop(valor)
        #print "Fitness pulento " +str(fm)
        #print "A Eliminar " + str(valor)



    return mejores

def seleccionespecial(population,umbral,datos):
    seleccionados = []

    for i in population:
        valor = 0
        for j in range (len(i.kmers)):
            for k in range (j,len(i.kmers)):
                distancia = HammingDistance(i.kmers[j].adn,i.kmers[k].adn)
                if (distancia==0 and j!=k):
                    valor = valor + 1
        if (valor>1):

            fitness= HammingMatrixr(i,datos)
            propor = decimal.Decimal(1/decimal.Decimal(fitness))

            i.setFitness(propor)
            seleccionados.append(i)

    print ("Seleccion Especial " +str(len(seleccionados)))

    sel = seleccion(seleccionados,umbral)
    return sel

def mutacion(individuo,datos):
    #print("MUTANDO")
    mutado = copy.deepcopy(individuo)
    cambios = random.randint(1,len(individuo.kmers)-1)
    ncambios = random.randint(1,len(individuo.kmers)-1)


    if (individuo.kmers[ncambios].posicioninicial + individuo.kmers[ncambios].largo + 1 > len(datos[ncambios])):
        kmer = Kmer(ncambios,individuo.kmers[ncambios].adn,individuo.kmers[ncambios].posicioninicial,individuo.kmers[ncambios].largo)
        kmer.posicioninicial = kmer.posicioninicial -1
        kmer.adn = getPalabra(datos,ncambios,kmer.posicioninicial,kmer.largo)
        individuo.setKmer(ncambios,kmer)
    else:
        kmer = Kmer(ncambios,individuo.kmers[ncambios].adn,individuo.kmers[ncambios].posicioninicial,individuo.kmers[ncambios].largo)
        kmer.posicioninicial = kmer.posicioninicial + 1
        kmer.adn = getPalabra(datos,ncambios,kmer.posicioninicial,kmer.largo)
        individuo.setKmer(ncambios,kmer)




    #for i in range(len(mutado.kmers)):
    #    start = random.randint(0,len(mutado.kmers)-1)
    #    kmer = Kmer(start,mutado.kmers[start].adn,mutado.kmers[start].posicioninicial,mutado.kmers[start].largo)
    #    posicion = random.randint(0,len(kmer.adn))
    #    kmer.posicioninicial = posicion
    #    kmer.adn = getPalabra(datos,start,kmer.posicioninicial,kmer.largo)
    #    mutado.setKmer(start,kmer)


    return individuo



def escribeGeneracion(population,archivo):
    f=open(archivo,'a')
    f.write("Population")
    f.write('\n')
    for i in population:
        valores = ""
        for j in i.kmers:
            valores = valores +"-"+ str(j.posicioninicial)
        f.write(valores)
        f.write('\n')
        f.write(str(i.fitness))
        f.write('\n')
    f.write('\n')
    f.write("---------")
    f.write('\n')
    f.close

def genetico(poblacion,ciclos,datos,lkmer):
    population = []
    resultados = []
    contador = 0
    fitnessacum = 0



    #Genera poblacion inicial
    for i in range (poblacion):
        m = generaMatrizInicial(datos,lkmer)
        #ftsd, pal = getFitness(m)


        population.append(m)


    pop = copy.deepcopy(population)

    for c in range(ciclos):
        promedio = 0
        suma = 0


        semuta = random.randint(0,100)
        #sin mutacion
        if (semuta<15):
            mutados = random.randint(0,len(pop))
            for m in range(mutados):

                amutar = random.randint(0,len(pop)-1)

                pop[amutar] = mutacion(pop[amutar],datos)
                fts,pal = getFitness(pop[amutar])
                propor = decimal.Decimal(decimal.Decimal(fts)/(decimal.Decimal(lkmer*len(datos))))
                pop[amutar].setFitness(propor)



        hijos = []
        poptemp = []
        for i in range(poblacion):
            cruza1 = random.randint(0,len(pop)-1)
            cruza2 = random.randint(0,len(pop)-1)
            #print "Cruza 1 " +str(cruza1)
            #print "Cruza 2 " +str(cruza2)
            hijo1,hijo2 = cruzamiento(pop[cruza1],pop[cruza2])

            pop.append(hijo1)
            pop.append(hijo2)

        #escribeGeneracion(pop,"Gens"+str(c)+".fts")
        #pop = hijos

        #pop = []
        #pop = copy.deepcopy(poptemp)





        q = Queue()
        tpop = len(pop)
        m1 = pop[:int(tpop/4)]
        p1 = Process(target=parallelFitness,args=(m1,q))
        p1.start()
        m2 = pop[int((tpop/4)+1):int((tpop/4)*2)]
        p2 = Process(target=parallelFitness,args=(m2,q))
        p2.start()
        m3 = pop[int(((tpop/4)*2)+1):int((tpop/4)*3)]
        p3 = Process(target=parallelFitness,args=(m3,q))
        p3.start()
        m4 = pop[int(((tpop/4)*3)+1):]
        p4 = Process(target=parallelFitness,args=(m4,q))
        p4.start()

        mejorfitness = 0
        presults = []
        pop = []
        mejorg = None
        for ps in range(4):
            presults.append(q.get(True))
        for pr in presults:
            for m in pr:
                if (m.fitness>mejorfitness):

                    mejorfitness = m.fitness
                    mejorg = copy.deepcopy(m)

                suma = suma + m.fitness
                pop.append(m)



        if (mejorg.fitness>fitnessacum):
            fitnessacum = mejorg.fitness




        promedio = suma/len(pop)

        largoantessel = len(pop)


        pop2 = copy.deepcopy(pop)
        qs = Queue()
        tpop = len(pop2)
        m1 = pop2[:int(tpop/4)]
        p1 = Process(target=parallelSeleccion,args=(m1,qs))
        p1.start()
        m2 = pop2[int((tpop/4)+1):int((tpop/4)*2)]
        p2 = Process(target=parallelSeleccion,args=(m2,qs))
        p2.start()
        m3 = pop2[int(((tpop/4)*2)+1):int((tpop/4)*3)]
        p3 = Process(target=parallelSeleccion,args=(m3,qs))
        p3.start()
        m4 = pop2[int(((tpop/4)*3)+1):]
        p4 = Process(target=parallelSeleccion,args=(m4,qs))
        p4.start()
        pop = []
        psseleccion = []

        for ps in range(4):
            psseleccion.append(qs.get(True))
        for pr in psseleccion:
            for m in pr:
                pop.append(m)


        #pop = seleccion(pop2,poblacion)
        #pop = seleccion(pop,umbral)

        escribeindividuo(mejorg,c,"results/"+sys.argv[1]+str(lkmer)+"-"+str(len(datos))+salida+"-generacion.fts")

        escribeFitness(promedio,lkmer,"results/"+sys.argv[1]+str(lkmer)+"-"+str(len(datos))+salida+"-fitnesspromedio.fts")

        print ("Generación " +str(c) + " Fitness Promedio " + str(promedio)+"Tamaño antes de sel "+str(largoantessel)+" Tamaño Poblacion "+str(len(pop))+" Mejor Individuop "+str(mejorg.fitness))


    #escribeGeneracion(pop,"GAFinal"+salida+".fts")

datos = cargaSecuencias(sys.argv[1])
salida = sys.argv[4]
#1 archivo, 2 ciclo
#poblacion,ciclos,datos,lkmer,umbral,salida,parada

    #annealing2(int(sys.argv[2]),int(sys.argv[3]),i,datos,salida)
lkmer = int(sys.argv[5])
genetico(int(sys.argv[2]),int(sys.argv[3]),datos,lkmer)
#for i in range(int(sys.argv[4]),int(sys.argv[5])):
        #ejecutar(int(sys.argv[2]),int(sys.argv[3]),i,datos,sys.argv[6],sys.argv[7])
        #annealing2(int(sys.argv[2]),int(sys.argv[3]),i,datos,sys.argv[6]+"a",sys.argv[7])
