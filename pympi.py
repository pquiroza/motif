import random
import copy
import decimal
import sys
import uuid
import hashlib
import time

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

archivo = "PS00010.fa"
poblacion = 25000
largo = 6
datos = []
pop = []


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


def showMatriz(matriz):
    for i in matriz.kmers:
        print (i.adn)
    print (matriz.fitness)




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



def parallelFitness(mmatriz):
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




if (rank==0):
    popglobal = []
    print("Soy el jefe")
    for i in range (1,size):
        lista = comm.recv(source=i)
        for l in lista:
            popglobal.append(l)

    for p in range(len(popglobal)):
        print(p,popglobal[p].fitness)


if (rank!=0):
    poplocal = []
    datos = cargaSecuencias(archivo)
    print(len(datos))
    for i in range(poblacion):
        m = generaMatrizInicial(datos,largo)
        poplocal.append(m)
    parallelFitness(poplocal)




    comm.send(poplocal,dest=0)
