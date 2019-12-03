import numpy as np

nodes = dict()  #Dictionary (Table in Excel)
elements = dict()

class Node:   #Object (Row in Excel)
  def __init__(self, ID, X, Y):  #constructor
    self.ID = ID
    self.X = X
    self.Y = Y
    self.rest = [0, 0, 0]    #restraints
    self.dof = [-1, -1, -1]  #dof number
    nodes[ID] = self


class Element:
  def __init__(self, ID, n1, n2, EI, EA):
    self.ID = ID
    self.n1 = nodes[n1]
    self.n2 = nodes[n2]
    self.EI = EI
    self.EA = EA


    Lx = self.n2.X - self.n1.X
    Ly = self.n2.Y - self.n1.Y
    L = (Lx**2 + Ly**2)**0.5
    c, s = Lx / L, Ly / L

    T = np.asarray(  
         [[c,  s, 0, 0, 0, 0],
         [-s, c, 0, 0, 0, 0],
         [0,  0, 1, 0, 0, 0],
         [0,  0, 0 ,c, s, 0],
         [0, 0, 0 ,-s, c, 0],
         [0, 0, 0,  0, 0 ,1]]  )


    KL = np.asarray(
      [[EA/L, 0, 0, -EA/L, 0, 0],
       [0, 12*EI/L**3, 6*EI/L**2, 0, -12*EI/L**3, 6*EI/L**2],
       [0, 6*EI/L**2, 4*EI/L, 0, -6*EI/L**2, 2*EI/L],
       [-EA/L, 0, 0, +EA/L, 0, 0],
       [0, -12*EI/L**3, -6*EI/L**2, 0, 12*EI/L**3, -6*EI/L**2],
       [0, 6*EI/L**2, 2*EI/L, 0, -6*EI/L**2, 4*EI/L]]
    ) 

    self.KG = T.T @ KL @ T
         
    elements[ID] = self

  def code(self):
    return  self.n1.dof + self.n2.dof


Node(ID=1, X=0, Y=0)   #Create (contruct) a node
Node(ID=2, X=2, Y=0)   #Create (contruct) a node
Node(ID=3, X=4, Y=0)   #Create (contruct) a node
Node(ID=4, X=6, Y=0)   #Create (contruct) a node
Node(ID=5, X=6, Y=-3)   #Create (contruct) a node

nodes[1].rest = [1, 1, 1]
nodes[5].rest = [1, 1, 0]

EI = 0.3*0.5**3/12*28e6
EA = 0.3 * 0.5 * 28e6

Element(ID=1, n1=1, n2=2, EI=EI, EA=EA)
Element(ID=2, n1=2, n2=3, EI=EI, EA=EA)
Element(ID=3, n1=3, n2=4, EI=EI, EA=EA)
Element(ID=4, n1=4, n2=5, EI=EI, EA=EA)

M = 0
for ID, n in nodes.items():
  if n.rest[0] == 0:
    n.dof[0] = M
    M = M + 1
  if n.rest[1] == 0:
    n.dof[1] = M
    M = M + 1
  if n.rest[2] == 0:
    n.dof[2] = M
    M = M + 1

N = M

for ID, n in nodes.items():
  if n.rest[0] == 1:
    n.dof[0] = M
    M = M + 1
  if n.rest[1] == 1:
    n.dof[1] = M
    M = M + 1
  if n.rest[2] == 1:
    n.dof[2] = M
    M = M + 1

print("N:", N, "M:", M)

for ID, n in nodes.items():
  print(ID,  n.X, n.Y, n.rest, n.dof)

KSistem = np.zeros((M, M))

for ID, e in elements.items():
  code = e.code()
  KElement = e.KG
  for i in range(6):
    for j in range(6):
      KSistem[code[i], code[j]] += KElement[i, j]



u1 = np.linalg.inv(KSistem[0:N, 0:N]) @ [0, -5, 0, 0, 0, -10, 0, 0 ,0 ,0]

print(u1)





  





