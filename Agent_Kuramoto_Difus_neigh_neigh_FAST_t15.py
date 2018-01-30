
# coding: utf-8

# In[25]:

#from pylab import *
import math
import random as rd
import numpy as np
import networkx as nx
from scipy import stats as st
from numpy import array
from scipy.integrate import ode
from scipy.stats import pareto
import csv

def sistema_nodo(nodo):  # genera las ecuaciones (derivada de cada P)
    
    global g
    
    resultado = []
    
    aux = str("g.node["+str(nodo)+"]['omega']")
    
    for i in g.neighbors_iter(nodo):
        
        aux += str("+alpha*np.sin(g.node["+str(i)+"]['theta']-P[0])") ###  P[]= fase (variable a integrar)
                
    resultado.append(aux)           
    
    return resultado

def dP_dt2_nodo(t,P,alpha,expresion,nodo):
    
    resultado = []
    
    resultado.append(eval(expresion[0]))           
    
    return resultado

def actualizacion2_nodo(g,Tg,nodo):
    
    global alpha
    
    expresion = sistema_nodo(nodo)
    #sistema
    #print g.nodes()
    #print(expresion)

    P0=[]

    P0.append(g.node[nodo]['theta'])

    t0 = 0.0
    t1 = Tg

    solver = ode(dP_dt2_nodo)
    solver.set_initial_value(P0, t0)
    solver.set_integrator('vode')
    solver.set_f_params(alpha,expresion,nodo)

    while solver.successful() and solver.t < t1:
        solver.integrate(t1, step=True)
   
    g.node[nodo]['theta'] = list(array(solver.y))[0]

def sistema_nodo_vecinos(nodo):  # genera las ecuaciones (derivada de cada P)
    
    global g
    
    candidatos = [nodo] + g.neighbors(nodo)
    
    #for index, s in enumerate(stocks_list):
    #print index, s
    
    resultado = []
    
    for i in candidatos: #la variable P de cada candidato será la de su índice en la lista: candidato[0]-->P[0]
    
        aux = str("g.node["+str(i)+"]['omega']")

        for j in g.neighbors_iter(i):
            
            if (j in candidatos):

                aux += str("+alpha*np.sin(P["+str(candidatos.index(j))+"]-P["+str(candidatos.index(i))+"])") ###  P[]= fase (variable a integrar)
                
            else:
                
                aux += str("+alpha*np.sin(g.node["+str(j)+"]['theta']-P["+str(candidatos.index(i))+"])") ###  P[]= fase (variable a integrar)

        resultado.append(aux)           
    
    return resultado

def dP_dt2_nodo_vecinos(t,P,alpha,expresion,nodo):
    
    resultado = []
    
    candidatos = [nodo] + g.neighbors(nodo)
    
    for i in range(0,len(candidatos)):
                
        resultado.append(eval(expresion[i]))            
    
    return resultado

def actualizacion2_nodo_vecinos(g,Tg,nodo):
    
    global alpha
    
    candidatos = [nodo] + g.neighbors(nodo)
    
    expresion = sistema_nodo_vecinos(nodo)
    
    #sistema
    #print g.nodes()
    #print(expresion)

    P0=[]
    
    for i in candidatos:

        P0.append(g.node[i]['theta'])

    t0 = 0.0
    t1 = Tg

    solver = ode(dP_dt2_nodo_vecinos)
    solver.set_initial_value(P0, t0)
    solver.set_integrator('vode')
    solver.set_f_params(alpha,expresion,nodo)

    while solver.successful() and solver.t < t1:
        solver.integrate(t1, step=True)
   
    for i in candidatos:

        g.node[i]['theta'] = list(array(solver.y))[candidatos.index(i)]

#########
def nuevo_nodo2(i): #Se selecciona un nodo y se le conecta a un nodo nuevo. Se perturba la posición del primero.
    
    global g, delta_omega, c, delta_c
     
    n_nodes = len(g.nodes())
    
    nuevo = int(n_nodes)
    g.add_node(nuevo)
    g.node[nuevo]['omega'] = 1.0 + np.random.uniform(-delta_omega, delta_omega)

    if n_nodes > 0:
  
        #i=rd.choice(range(n_nodes)) #J: se selecciona el nodo en el update2()
        
        g.add_edge(nuevo, i)
        
        g.node[nuevo]['theta'] = g.node[i]['theta'] + np.random.uniform(-delta_theta, delta_theta) #A: nace con una fase similar a la del
        #"padre" (se introduce una perturbacion)
        angulo = rd.random()*2*math.pi
        modulo = rd.random()*c
        g.node[nuevo]['x']=g.node[i]['x'] + modulo*np.cos(angulo)
        g.node[nuevo]['y']=g.node[i]['y'] + modulo*np.sin(angulo)
        
        #angulo = rd.random()*2*math.pi
        #modulo = rd.random()*delta_c
        #g.node[i]['x'] += modulo*np.cos(angulo)
        #g.node[i]['y'] += modulo*np.sin(angulo)
        
        return i #devuelve el valor del nodo al que se ha conectado el nodo recien creado
    
    else:
        
        g.node[nuevo]['theta'] = 2.0 *math.pi * rd.random()
        g.node[nuevo]['x'] = rd.random()
        g.node[nuevo]['y'] = rd.random()
        
        return nuevo #devuelve el valor del nodo al que se ha conectado el nodo recien creado


def nuevo_enlace2(i): #Se selecciona un nodo y se le conecta a otro cercano. Se perturba la posición del primero.
    
    global g, delta_omega, b, delta_c
    
    max_intentos = 1000
    exito = False
    
    n_nodes = len(g.nodes())
    n_enlaces = len(g.edges())
    max_enalces = 0.5 * n_nodes * (n_nodes - 1.0)
    x0=g.node[i]['x']
    y0=g.node[i]['y']    
    
    next_neighbors=[] #J:  primero busco los 2º vecinos
    for ii in nx.neighbors(g,i):
        for nnb in nx.neighbors(g,ii):
            #print('i=',i,'ii=',ii)
            next_neighbors.append(nnb)  #J: Solo consideramos enlaces con los 2º vecinos
            #print('next_neigh=',next_neighbors)
            
    
    if (n_enlaces < max_enalces) and (len(next_neighbors)>0):
        
        #neighbs = nx.neighbors(g,i)
        intentos = 0 
        while (exito != True) and (intentos<max_intentos):
            intentos += 1 #J: Hay que dejar siempre una vía de escape en los while
            #i = rd.choice(g.nodes()) #El agente se elige en el update2()

            rad = pareto.rvs(b, size=1) #Generate random numbers form a pareto density distribution b/(x^(1+b))
            rad = rad[0] #El comando anterior genera una lista con 1 elemento, que extraemos aquí

            #candidates0 = [nb for nb in g.nodes() if ( ((g.node[nb]['x']-x0)**2 + (g.node[nb]['y']-y0)**2) < rad**2) and (nb != i) ]
            candidates0 = [nb for nb in next_neighbors if ( ((g.node[nb]['x']-x0)**2 + (g.node[nb]['y']-y0)**2) < rad**2) and (nb != i) ]
            
            candidates = [nb for nb in candidates0 if (nb in nx.non_neighbors(g,i))] #J: antes de elegir vecino comprueba que no tiene enlace
            n_candidates = len(candidates)  #se limita por los de su especie y la otra
            #print('i:',i,'n_candidates= ',n_candidates)

            if n_candidates > 0:

                j=rd.choice(candidates)

                if j in nx.non_neighbors(g,i): #nx.non_neighbors(g,i) proporciona la lista de no-vecinos de i
                    g.add_edge(i, j)
                    exito = True
                    
            #angulo = rd.random()*2*math.pi
            #modulo = rd.random()*delta_c
            #g.node[i]['x'] += modulo*np.cos(angulo)
            #g.node[i]['y'] += modulo*np.sin(angulo)
            
        return i #devuelve el valor del nodo al que se ha conectado el nodo recien creado
                    
    else: #Si el grafo es completo, añadimos un nuevo nodo
        
        #i = rd.choice(g.nodes()) #J: Si no puede añadir enlace VUELVE
        #nuevo_nodo2(i)
        return i
    


def update2():
    
    global g, g_max, n_max, t, umbral, Tg, m, r, delta_omega, p, temp, p_rotura, gamma, exp_Levy,delta_c,delta_theta
    
        
    ##################################################    #Adición de un nuevo nodo o enlace (según sea el valor de p)  
    nodo = rd.choice(g.nodes()) #J: Se elige primero el agente (Sacamos un candidatos de la lista)

    angulo = rd.random()*2*math.pi #J: 1º Hacemos difusión del nodo con amplitud Delta_c
    modulo = rd.random()*delta_c
    g.node[nodo]['x'] += modulo*np.cos(angulo)
    g.node[nodo]['y'] += modulo*np.sin(angulo)
    
    
    
    if rd.random()<= p:
        
        #nodo = nuevo_enlace2(i) 
        nuevo_enlace2(nodo) #J: Ahora de devuelve nada, yo le mando el nodo

    else:
        
        #nodo = nuevo_nodo2(i)
        nuevo_nodo2(nodo) #J: Ahora de devuelve nada, yo le mando el nodo

##################################################    #Se da Tg microtiempos para np.sincronizarse
    if(incluye_vecinos):
        actualizacion2_nodo_vecinos(g,Tg,nodo) #En cada paso sólo se modifican las fases del nodo donde se añade el nuevo nodo/enlace
    else:
        actualizacion2_nodo(g,Tg,nodo)
    #y la de sus primeros vecinos
    
######################################################  #Evaluación de la estabilidad / Relajación     

    escision = 0 #Variable que almacena el tamaño de nodos escindidos en un instante de tiempo
    rotura = 0
            
    revisar = True
    
    while revisar == True:
    
        lista=[] #Almaceno inestables

        if len(g.edges()) > 0:
             
        
            for i in g.neighbors_iter(nodo):
                
                media_sin = 0.5 * (np.sin(g.node[nodo]['theta']) + np.sin(g.node[i]['theta']))
                media_cos = 0.5 * (np.cos(g.node[nodo]['theta']) + np.cos(g.node[i]['theta']))
                
                r_ij = abs(complex(media_cos,media_sin)) #Parámetro de orden del enlace
                
                d2_ij= (g.node[nodo]['x']-g.node[i]['x'])**2 + (g.node[nodo]['y']-g.node[i]['y'])**2
                
                d_ij = math.sqrt(d2_ij) #J: distancia=r NO r^2
                
                
                if d_ij > 1:
                    
                    umbral_sincro_d = 1.0 - (1.0 - umbral)/(d_ij**(gamma)) 
                    p_rotura_d = 1.0 - (1.0 - p_rotura)/(d_ij**(gamma))
                    
                else:
                    umbral_sincro_d = umbral 
                    p_rotura_d = p_rotura
                
                #print('(i:',i[0],',j:',i[1],');d_ij1000 = {:.4f}'.format(d_ij1000),'umbral={:.4f}'.format(umbral_sincro_d), ' r_ij={:.4f}'.format(r_ij),' p_rotura={:.4f}'.format(p_rotura_d))

                if  (r_ij < umbral_sincro_d) or (rd.random() <= p_rotura_d):

                    lista.append(i)
                    #print('lista de enlaces a cortar',lista)

            if lista == []:

                revisar = False

            else:

                for i in lista:
                    
                    g.remove_edge(nodo,i) #Eliminamos el enlace
                    
                    rotura +=1
                    
                    #for j in [i[0],i[1]]: #Si tras la rotura algún nodo se queda desconectado, Levy-flight
                        
                    #    if g.degree(j) == 0:
                            
                    #        #print('nodo desconectado:',j)
                    #        angulo = rd.random()*2*math.pi
                    #        salto = (rd.random())**(-1/exp_Levy)

                    #        g.node[j]['x'] += salto*np.cos(angulo)
                    #        g.node[j]['y'] += salto*np.sin(angulo)
                    #        #print(g.node[j])
                            
                            
                
        else:
            
            revisar = False
            
    lista_roturas.append(rotura)
    
    
def registro_datos():
    
    global g, lista_cant_nodos,lista_cant_enlaces,lista_roturas,lista_media_sin,lista_media_cos,    lista_orden,lista_cant_subgrafos,lista_cantidad_nodos_componente_gigante,lista_diametro_componente_gigante,    lista_clustering_medio_componente_gigante,lista_camino_corto_medio_componente_gigante,    lista_distancia_enlace_euclidiana_media_componente_gigante
    
    lista_cant_nodos.append(len(g.nodes()))
    lista_cant_enlaces.append(len(g.edges()))
    #lista_escisiones.append(escision)
    #lista_roturas.append(rotura) #A: las roturas se registran en la funcion update,\
    #para así tener el dato de los instantes de microtiempo
    
    diferencias_sin = 0.0
    diferencias_cos = 0.0
        
    for i in g.nodes_iter():
            
        diferencias_sin += np.sin(g.node[i]['theta'])/len(g.nodes())
        diferencias_cos += np.cos(g.node[i]['theta'])/len(g.nodes())
 
    lista_media_sin.append(diferencias_sin)
    lista_media_cos.append(diferencias_cos)
    lista_orden.append(abs(complex(diferencias_cos,diferencias_sin)))
    
    ##################################################    #Registro datos sobre los subgrafos y la componente gigante
    
    subgrafos = sorted(nx.connected_components(g), key = len, reverse=True)
    
    lista_cant_subgrafos.append(len(subgrafos))
    
    componente_gigante = g.subgraph(list(subgrafos[0])) #A: Grafo que contiene una copia del componente gigante
    
    lista_cantidad_nodos_componente_gigante.append(len(componente_gigante.nodes())) #A: Almacena la variable que indica el nombre
    lista_diametro_componente_gigante.append(nx.diameter(componente_gigante)) #A: Almacena la variable que indica el nombre
    
    if len(componente_gigante.nodes()) > 1:
        
        lista_clustering_medio_componente_gigante.append(nx.average_clustering(componente_gigante,weight=None))
        lista_camino_corto_medio_componente_gigante.append(nx.average_shortest_path_length(componente_gigante, weight=None))
        
    else:
        
        lista_clustering_medio_componente_gigante.append(float('nan'))
        lista_camino_corto_medio_componente_gigante.append(float('nan'))
        
    if len(componente_gigante.edges()) > 1:
        
        distancias_aux = []
        
        for i in componente_gigante.edges():
            
            d2_ij= (g.node[i[0]]['x']-g.node[i[1]]['x'])**2 + (g.node[i[0]]['y']-g.node[i[1]]['y'])**2
            
            distancias_aux.append(math.sqrt(d2_ij))
    
        lista_distancia_enlace_euclidiana_media_componente_gigante.append(1.0*sum(distancias_aux)/len(distancias_aux))
    
    else:
        
        lista_distancia_enlace_euclidiana_media_componente_gigante.append(float('nan'))
    
#############################################################
#############################################################
#PARAMETROS
#############################################################
#Parámetros generales
tmax = 15
p = 0.5 # parametro para ajustar la proporcion de nodos o enlaces nuevos que introducimos en cada paso
incluye_vecinos = True #J: Parámetro para elegir si actualizar también las fases de los vecinos o no

#a = 0.001; # limiting factor
b = 2.0  #Pareto parameter
c = 2.0 #factor distancia nuevo hijo  (por ahora no se mueren)
gamma = 2.0 #exponente de Pareto
exp_Levy = 1.5 #exponente de Levy-Flights

delta_omega = 0.2
delta_theta = 0.1*math.pi #máxima desviación de la fase del hijo
delta_c = 0.1*c ##máxima desviación de la posición de un nodo al inicio de un nuevo instante de tiempo

#############################################################
#Umbrales

umbral = 0.75 #sincronización mínima que se exige en los extremos de un enlace
p_rotura = 0.0 #probabilidad de rotura de un enlace que está a distancia menor o igual que 1

#############################################################
#Integración numérica
Tg = 3.0 #duración de la integración numérica de las ecuaciones diferenciales
alpha = 1.0 # coupling strength

csvfile = 'salida_listas_Agentes_nn.txt' #nombre del fichero txt de salida con los datos de las series temporales
path_graph = 'salida_grafo_Agentes_nn.gpickle' #nombre del fichero gpickle de salida con el grafo

lista_cant_nodos = []
lista_cant_enlaces = []
#lista_escisiones = []
lista_roturas = []
lista_media_sin = []
lista_media_cos = []
lista_orden = []
lista_cant_subgrafos = []
lista_cantidad_nodos_componente_gigante = [] #A: Almacena la variable que indica el nombre
lista_diametro_componente_gigante = [] #A: Almacena la variable que indica el nombre
tiempo = []
microtiempo = [] #Almacena el valor de (macro)tiempo de cada instante de microtiempo
lista_clustering_medio_componente_gigante = [] #A: Almacena la variable que indica el nombre
lista_camino_corto_medio_componente_gigante = [] #A: Almacena la variable que indica el nombre
lista_distancia_enlace_euclidiana_media_componente_gigante = []  #A: Almacena la variable que indica el nombre

##################################################  Semilla

g = nx.Graph()
nuevo_nodo2(0)

#g_max = g.copy()
#n_max = len(g_max.nodes())

t = 1

lista_cant_nodos.append(1)
lista_cant_enlaces.append(0.0)
#lista_escisiones.append(0.0)
lista_roturas.append(0.0)
lista_cant_subgrafos.append(1.0)
lista_media_sin.append(np.sin(g.node[0]['theta']))
lista_media_cos.append(np.cos(g.node[0]['theta']))
lista_orden.append(1.0)
lista_cantidad_nodos_componente_gigante.append(1.0)
lista_diametro_componente_gigante.append(float('nan'))
tiempo.append(0)
microtiempo.append(0)
lista_clustering_medio_componente_gigante.append(float('nan'))
lista_camino_corto_medio_componente_gigante.append(float('nan'))
lista_distancia_enlace_euclidiana_media_componente_gigante.append(float('nan'))

################################
################################
## Comienzo del bucle
################################

cons = 0

while t < tmax+1:
    
    dt=0
    cons = 0
    
    while dt < 1.0:
        
        cons += 1
        
        if cons % 1000 == 0:
        
            print('cons=',cons)
        
        microtiempo.append(t)
        
        dt += 1./len(g.nodes())
        #print('t:',t,';  n_nodes=',len(g.nodes()),'; dt: ',dt)
        update2() #J: A continuación habría que meter lo que se calcula en cada actualización
        
    registro_datos() #Esta función registra los datos en los instante de macro tiempo
    tiempo.append(t)
    t += 1

    #J: Uso la actualización de tiempos de SAYAMA
    #J: Aquí habría que meter lo que se calcula para cada MACROTIEMPO
    
    if t%1==0:
        
        print(t)
        res = [lista_cant_nodos,lista_cant_enlaces,lista_roturas,lista_media_sin,lista_media_cos,lista_orden,               lista_cant_subgrafos,lista_cantidad_nodos_componente_gigante,lista_diametro_componente_gigante,tiempo,               microtiempo,lista_clustering_medio_componente_gigante,lista_camino_corto_medio_componente_gigante,               lista_distancia_enlace_euclidiana_media_componente_gigante]

        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(res)
        
        nx.write_gpickle(g, path_graph)

#Salvado final, tras completar la ejecución del programa

res = [lista_cant_nodos,lista_cant_enlaces,lista_roturas,lista_media_sin,lista_media_cos,lista_orden,       lista_cant_subgrafos,lista_cantidad_nodos_componente_gigante,lista_diametro_componente_gigante,tiempo,       microtiempo,lista_clustering_medio_componente_gigante,lista_camino_corto_medio_componente_gigante,       lista_distancia_enlace_euclidiana_media_componente_gigante]

with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(res)
    
nx.write_gpickle(g, path_graph) 


# In[21]:

print(nx.neighbors(g,4))

