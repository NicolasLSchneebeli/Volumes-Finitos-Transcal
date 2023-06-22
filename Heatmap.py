import matplotlib.pyplot as plt
import numpy as np

#Código feito por Arthur Frasson Pretti e Nicolas Lovatte Schneebeli em 11/11/2022 para avaliação em Transferência de Calor I em UFES 2022/2

#Dimensões
lx1=5e-3
ly1=3e-3
ly0= 2e-3 #tamanho da pá em cima da parte com Ar
#Carregando o arquivo
temp=np.loadtxt("tempfinal15vol.txt")
tamanho_temp=temp.shape

x, y = np.meshgrid(np.linspace(0,lx1,tamanho_temp[1]), np.linspace(ly1,0,tamanho_temp[0]))

z = temp.copy()
z=np.flip(z,0)



zmin=np.abs(np.nanmin(z))   
zmax=np.abs(np.nanmax(z))

fig, ax = plt.subplots()
lz=np.logical_not(np.isnan(z))
mz=np.ma.masked_where(np.isnan(z),z)
print(mz)

c = ax.pcolormesh(x, y, mz, cmap='hot_r',vmin=zmin, vmax=zmax) #,shading='gouraud')
ax.set_title('Mapa de calor')
plt.xlabel("x[m]")
plt.ylabel("y [m]")
ax.axis([x.min(), x.max(), y.min(), y.max()])
plt.gca().invert_yaxis()
fig.colorbar(c, ax=ax)
plt.show()

#Fazendo um corte entre os gases de combustão e o ar
corte=temp[:,int(3*tamanho_temp[1]/5)]
corte=corte[np.logical_not(np.isnan(corte))]

y1=np.linspace(0,ly0,int(2*tamanho_temp[0]/3+1))
plt.plot(y1,corte,"o",ls='-')
plt.xlabel("Y[m]")
plt.ylabel("Temperatura[K]")
plt.show()

#Fluxo 

u=np.loadtxt("fluxofinalu15vol.txt") #Abrindo arquivos da simulação
v=np.loadtxt("fluxofinalv15vol.txt")
u1=u.round(decimals=2) #Arredondamento para melhor visualização
v1=v.round(decimals=2)


#Gerando campo vetorial:
x2=np.linspace(0,ly1,tamanho_temp[0])
y2=np.linspace(0,lx1,tamanho_temp[1])
x2,y2=np.meshgrid(y2,x2)
plt.quiver(x2,y2,u,-v)
plt.gca().invert_yaxis()
plt.show()


c = ax.pcolormesh(x, y, mz, cmap='hot_r',vmin=zmin, vmax=zmax) #,shading='gouraud')
ax.axis([x.min(), x.max(), y.min(), y.max()])
plt.gca().invert_yaxis()
fig.colorbar(c, ax=ax)
plt.quiver(x2,y2,u,-v)
plt.xlabel("x[m]")
plt.ylabel("y [m]")
ax.set_title('Mapa de calor')
plt.show()