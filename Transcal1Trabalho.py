import matplotlib.pyplot as plt
import numpy as np
import time


#Código feito por Arthur Frasson Pretti e Nicolas Lovatte Schneebeli em 11/11/2022 para avaliação em Transferência de Calor I em UFES 2022/2

tempo = time.time() #Contar o tempo de execução da simualção

# Dimensoes
lx1=5e-3
lx0=2e-3
ly1=3e-3
ly0=2e-3

#Dados de criação
nvoly=15#volumes em x (Ideal é ser múltiplo de 15)
nvolx=int(lx1*nvoly/ly1) #volumes em y

# Regras de iteracao
itermax=8000000
tol=1.0e-6
tolmax=1.0e2/tol

# Volumes finitos tem que ser iguais
dy= ly1/nvoly
dx= lx1/nvolx

assert abs(dy-dx) < tol, "Dx diferente de Dy" #Se isso acontecer hipotese de dx=dy não é válida
    

#Dados do problema
cond=25 #W/mK
T_ar_int=400 #K
h_int=200 #W/m2K
T_gases_ext= 1700 #K
h_gas_ext=1000 #W/m2K

#Definindo volumes de controle vide relatório
ivertice = int(ly0*nvoly/ly1+0.5) #N° do Volume do vértice com convecção interna 
jvertice = int(lx0*nvolx/lx1+0.5)

assert abs(dx*jvertice-lx0) < tol, "Discretização em X ruim"
assert abs(dy*ivertice-ly0) < tol, "Discretização em Y ruim"



t=np.zeros([nvoly+1,nvolx+1])
m=np.copy(t) #Matriz que determinará as equações vigentes

#Preenchendo m:


m[0][1:] = 1    #Código: quando m=1 há condução em 3 pontos e convecção na Parte superior

m[:,ivertice][jvertice:]=2  #Código: quando m=2 há condução em 3 pontos e Convecção com o ar na Lateral Esquerda

m[ivertice][jvertice:]=3    #Código: quando m=3 há condução em 3 pontos e convecção na Parte Inferior
 
m[:,0][0:]=4    #Código: quando m=4 há 3 Conduções e 1 isolamento na Lateral Esquerda

m[:,nvolx][1:ivertice]=5  #Código: quando m=5 há Isolamento na Lateral Direita e 3 Conduções

m[nvoly][1:jvertice]=6  #Código: quando m=6 há Isolamento na Parte Inferior e 3 Conduções

m[m==0]=7   #Código: quando m=7 há 4 conduções

m[ivertice][jvertice]=9 #Código m: 9, ponto nodal em um vértice interno com convecção (Quina na pá)

m[:,ivertice+1:nvolx+1][jvertice+1:nvoly+1]=0   #Preenchendo o Ar

#Casos Especiais 
m[0][0]=8   #Canto Superior Esquerdo
m[0][nvolx]=8   #Canto Superior Direito
m[nvoly][0]=8   #Canto inferior direito
m[ivertice][nvolx]=8    #Canto direito onde o Ar intersecta a pá
m[nvoly][jvertice]=8    #Canto inferior onde o ar intersecta a par

iter=0
erro=1
hist_erro=[]

#Com as eq do incropera CAPÍTULO 4 TABELA 4.2 5ª Ed Adaptado.
while (iter <= itermax) and (erro >=tol) and (erro<=tolmax):
    told=np.copy(t)
    for i in range (0,nvoly+1):
        for j in range (0,nvolx+1):
            if m[i][j]==0:
              continue  
            elif m[i][j]==1: 
                t[i][j]=(2*t[i+1][j]+t[i][j-1]+t[i][j+1]+2*h_gas_ext*T_gases_ext*dy/cond)/(2*h_gas_ext*dy/cond+4)
                
            elif m[i][j]==2:
                 t[i][j]=(2*t[i][j-1]+t[i+1][j]+t[i-1][j]+2*h_int*T_ar_int*dx/cond)/(2*h_int*dx/cond+4)   
                 
            elif m[i][j]==3:
                t[i][j]=(2*t[i-1][j]+t[i][j-1]+t[i][j+1]+2*h_int*T_ar_int*dx/cond)/(2*h_int*dx/cond+4)   
                
            elif m[i][j]==4:
                t[i][j]=(t[i-1][j]+t[i+1][j]+2*t[i][j+1])/4
                
            elif m[i][j]==5:
                t[i][j]=(t[i+1][j]+t[i-1][j]+2*t[i][j-1])/4
                
            elif m[i][j]==6:
                t[i][j]=(t[i][j+1]+2*t[i-1][j]+t[i][j-1])/4
                
            elif m[i][j]==7:
                t[i][j]=(t[i-1][j]+t[i+1][j]+t[i][j-1]+t[i][j+1])/4
                
            elif m[i][j]== 8: #Médias dos pontos vizinhos: não afeta tanto a simulação. 
                
                t[0][0]=((h_gas_ext * dx/cond) * T_gases_ext + t[0][1] + t[1][0])/(h_gas_ext * dx/cond + 2) #OK
                
                t[0][nvolx]=((h_gas_ext * dx/cond) * T_gases_ext + t[0][nvolx-1] + t[1][nvolx])/(h_gas_ext * dx/cond + 2) #OK
                
                t[nvoly][0]=(t[nvoly-1][0]+t[nvoly][1])/2 #OK
                
                t[ivertice][nvolx]=((h_int*dx/cond)*T_ar_int+(t[ivertice][nvolx-1]+t[ivertice-1][nvolx]))/(h_int*dx/cond +2)
                
                t[nvoly][jvertice]=((h_int*dx/cond)*T_ar_int +(t[nvoly][jvertice-1]+t[nvoly-1][jvertice]))/(h_int*dx/cond +2)
                
            elif m[i][j]==9:
                t[i][j]=(2*(t[i][j-1]+t[i-1][j]) + (t[i][j+1]+t[i+1][j]) + 2*dx*h_int*T_ar_int/cond)/(6+2*h_int*dx/cond)                      
    diferenca=np.absolute(t-told)
    maior=np.max(diferenca)
    erro=maior*maior #Definição do Erro
    print("a iteração é {} e o erro {}".format(iter, erro))
    hist_erro.append(erro)
    iter=iter+1
    #Visualizando o motivo de parada!
    if iter >= itermax:
        print("Atingiu n° de iterações máx")
        
    if erro <= tol:
        print("Atingiu tolerância mín")
        
    if erro >= tolmax:
        print("Atingiu tolerância máx")
        
Temperatura_final=t.round(decimals=2) #Arredondamento para melhor visualização
print(Temperatura_final)

t[t==0]=np.nan #Retirando os 0, vai ser útil mais pra frente para calcular os mín e máx para o heatmap
np.savetxt("tempfinal{}vol.txt".format(nvoly),t) #Para salvar em um txt em outro execução faremos o heatmap
# #Favor abrir o arquivo .py chamado de "Heatmap.py" se quiser ver apenas do heatmap,campo vetorial, corte e dos erros sem a simulação 

#Plotando o erro pelo n° de iterações em log
iterlin=np.linspace(0,iter,iter)
plt.plot(np.log(iterlin),hist_erro)
plt.xlabel("N° de iterações(Log)")
plt.ylabel("Tamanho do erro")
plt.show()


#Fazendo um corte entre os gases de combustão e o ar
corte=t[0][:]
x1=np.linspace(0,lx1,nvolx+1)
plt.plot(x1,corte,"o",ls='-')
plt.xlabel("x[m]")
plt.ylabel("Temperatura[K]")
plt.show()

# #Cálculo do Fluxo

u=np.zeros([nvoly+1,nvolx+1])
v=np.copy(u)


# #Calculando o fluxo x e y:

for i in range (0, nvoly+1):
    for j in range (0,nvolx+1):
        
        if m[i][j] == 0:
            continue
        elif m[i][j]==1:
            u[i][j] =   -cond*(t[i][j + 1] - t[i][j]) / dx  #Fluxo na direção x
            v[i][j] =   -cond*(t[i+1][j]-t[i][j])/dy    #Fluxo na direção y

        elif m[i][j]==2:
            u[i][j] =   h_int*(t[i][j]-T_ar_int)    #Fluxo na direção x
            v[i][j] =   -cond * (t[i + 1][j] - t[i][j]) / dy    #Fluxo na direção y

        elif m[i][j]==3:
            u[i][j] =   -cond * (t[i][j + 1] - t[i][j]) / dx   #Fluxo na direção x
            v[i][j] =   h_int*(t[i][j]-T_ar_int)    #Fluxo na direção y
 

        elif m[i][j]==4:
            u[i][j] =   -cond * (t[i][j + 1]-t[i][j]) / dx  #Fluxo na direção x
            v[i][j] =   -cond * (t[i + 1][j] -t[i][j]) / dy     #Fluxo na direção y

        elif m[i][j]==5:
            u[i][j] =   0  ## Fluxo na direção x
            v[i][j] =   -cond*(t[i + 1][j] - t[i][j])/dy   #Fluxo na direção y

        elif m[i][j]==6:
            u[i][j] =   -cond*(t[i][j + 1] - t[i][j])/dx    #Fluxo na direção x
            v[i][j] =   0   #Fluxo na direção y

        elif m[i][j]==7:
            u[i][j]=    -cond*(t[i][j+1]-t[i][j])/dx    #Fluxo na direção y
            v[i][j]=    -cond*(t[i+1][j]-t[i][j])/dy    #Fluxo na direção y

        elif m[i][j]==8:
        #Ponto 0,0
            u[0][0]=    -cond*(t[0][1]-t[0][0])/dx
            v[0][0]=    -cond*(t[1][0]-t[0][0])/dy
        
        #Ponto 0, nvolx
            u[0][nvolx]=    0
            v[0][nvolx]=    -cond*(t[1][nvolx]-t[0][nvolx])/dy
        
        #Ponto nvoly,0
            u[nvoly][0] = -cond*(t[nvoly][1]-t[nvoly][0])/dx
            v[nvoly][0]=  0
        #Ponto 2/3 nvoly, nvolx
            u[ivertice][nvolx] =    0
            v[ivertice][nvolx] =    h_int*(t[ivertice][nvolx]-T_ar_int)
            
        #Ponto nvoly, 2/5 nvolx 
            u[nvoly][jvertice] =    h_int*(t[nvoly][jvertice]-T_ar_int)
            v[nvoly][jvertice] =    0   
        elif m[i][j]==9:
            u[i][j]=    h_int*(t[i][j]-T_ar_int)-cond*(t[i][j+1]-t[i][j])  
            v[i][j]=    h_int*(t[i][j]-T_ar_int)-cond*(t[i+1][j]-t[i][j])  

np.savetxt("fluxofinalu{}vol.txt".format(nvoly), u) #Salvando os arquivos de u e v
np.savetxt("fluxofinalv{}vol.txt".format(nvoly),v)


#Gerando campo vetorial e heatmap

#Heatmap
x1, y1 = np.meshgrid(np.linspace(0,lx1,nvolx+1), np.linspace(ly1,0,nvoly+1))
z = t.copy()
z=np.flip(z,0)

zmin=np.abs(np.nanmin(z))   
zmax=np.abs(np.nanmax(z))

fig, ax = plt.subplots()


c = ax.pcolormesh(x1, y1, z, cmap='hot_r',vmin=zmin, vmax=zmax)
ax.set_title('Mapa de calor')
plt.xlabel("x[m]")
plt.ylabel("y[m]")
ax.axis([x1.min(), x1.max(), y1.min(), y1.max()])
plt.gca().invert_yaxis()
fig.colorbar(c, ax=ax)
plt.show()


#Campo vetorial
x=np.linspace(ly1,0,nvoly+1)
y=np.linspace(0,lx1,nvolx+1)
x,y=np.meshgrid(y,x)
plt.quiver(x,y,u,v)
plt.gca().invert_yaxis()
plt.show()
print("--- {} segundos de simulação ---".format(time.time() - tempo)) #Tempo final de simulação
