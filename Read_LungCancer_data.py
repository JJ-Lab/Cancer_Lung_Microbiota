#!/usr/bin/env python
# coding: utf-8

## Software to read lung cancer data 
## data: 12/01/2020 
## authors: Juan Manuel Pastor and Javier Galeano

import pandas as pd
import networkx as nx
import numpy as np



# #### Reading data and metadata



df_data = pd.read_excel('Data/tabla-cancer-pulmon2.xlsx', header=1)
df_meta = pd.read_excel('Data/metadata-cancer-pulmon2.xlsx', header=0)

# Use df_data.head() or df_data.tail() to control the data


#### Divided dataset in Healthy and Sick samples 


df_h=df_meta[df_meta.loc[:,'Clinical']=='healthy']
df_s=df_meta[df_meta.loc[:,'Clinical']=='sick']

# df_s.head()


# Another way to divide the dataset 


df_s2=df_meta[df_meta['Clinical']=='sick']
df_h2=df_meta[df_meta['Clinical']=='healthy']

# df_s2.tail()


# ### We need taxa names to use as columns
# 


#taxa=df_data.loc[:,['#TAXA ID']]
#taxa0=taxa.squeeze()
#taxa2=[]
#for i in range(len(taxa0)):
#    family = taxa0[i].split('|')
#    taxa2.append(family[-1])
#    
#len(taxa2)

# To make different strangers cases as "uncultured", "Ambiguous",
#  "__" and so on, wa add a control number

taxa=df_data.loc[:,['#TAXA ID']]
taxa0=taxa.squeeze()
taxa1=[]
for i in range(len(taxa0)):
    family = taxa0[i].split('|')
    if (family[-1]=='__' ) or ( 'uncultured' in family[-1]) or ('Ambiguous' in family[-1]) or ('metagenome' in family[-1]) or ('unidentified' in family[-1]):
        #taxa1.append('?')
        taxa1.append(family[-2]+str(i))
    else:
        taxa1.append(family[-1])


# ### We divide our data in 6 different types:
# Healthy
# 1. saliva
# 2. healthy-lung
# Sick
# 1. saliva
# 2. healthy-lung
# 3. affected-lung
# 4. faeces



# Healthy data
df_h2_saliva=df_h2[df_h2['Location-longit']=='saliva']
df_h2_hlung=df_h2[df_h2['Location-longit']=='healthy-lung']

data_h_saliva=df_data[df_h2_saliva['sample']]
data_h_hlung=df_data[df_h2_hlung['sample']]
df_h_saliva =pd.DataFrame(data=data_h_saliva.values,index=taxa1,columns=data_h_saliva.columns)
df_h_hlung =pd.DataFrame(data=data_h_hlung.values,index=taxa1,columns=data_h_hlung.columns)



# Sick data
df_s2_saliva=df_s2[df_s2['Location-longit']=='saliva']
df_s2_hlung=df_s2[df_s2['Location-longit']=='healthy-lung']
df_s2_alung=df_s2[df_s2['Location']=='affected-lung']
df_s2_faeces=df_s2[df_s2['Location']=='faeces']


data_s_saliva=df_data[df_s2_saliva['sample']]
data_s_hlung=df_data[df_s2_hlung['sample']]
data_s_alung=df_data[df_s2_alung['sample']]
data_s_faeces=df_data[df_s2_faeces['sample']]

df_s_saliva =pd.DataFrame(data=data_s_saliva.values,index=taxa1,columns=data_s_saliva.columns)
df_s_hlung =pd.DataFrame(data=data_s_hlung.values,index=taxa1,columns=data_s_hlung.columns)
df_s_alung =pd.DataFrame(data=data_s_alung.values,index=taxa1,columns=data_s_alung.columns)
df_s_faeces =pd.DataFrame(data=data_s_faeces.values,index=taxa1,columns=data_s_faeces.columns)


# ### Save dataset as csv 
# 

df_s_faeces.to_csv('faeces_s.csv',index =True, header=True)
df_s_saliva.to_csv('saliva_s.csv',index =True, header=True)
df_s_hlung.to_csv('hlung_s.csv',index =True, header=True)
df_s_alung.to_csv('alung_s.csv',index =True, header=True)

df_h_saliva.to_csv('saliva_h.csv',index =True, header=True)
df_h_hlung.to_csv('hlung_h.csv',index =True, header=True)


### Building the network
### We define a function to build-up the networks.

def make_graph_median(graph):
    ma_graph=nx.Graph()
    Dic={}
    
    for patien in graph.columns:
        print(patien)
        graph_0 = graph[graph[patien]!=0][patien]
        suma = np.sum(graph_0)
        s16S_norm=graph_0/suma
        
        g_m = nx.complete_graph(len(graph_0))
        
        labels={}
        for i, tax_i in zip(range(len(graph_0)),graph_0.index):
            labels[i]=tax_i
            g_m=nx.relabel_nodes(g_m,labels)
            
        for u,v,d in g_m.edges(data=True):
            w_uv = s16S_norm[u]*s16S_norm[v]
                
            if (u,v) not in ma_graph.edges():
                ma_graph.add_edge(u,v,weight=[w_uv])
                Dic.setdefault((u,v),[])
                Dic[(u,v)].append(w_uv)
            
            else:
                Dic[(u,v)].append(w_uv)
                ma_graph[u][v].update(weight=Dic[(u,v)])

    for u,v,d in ma_graph.edges(data=True):
        dato=d['weight']
        ma_graph[u][v].update(weight=np.median(dato))


    return ma_graph
        
def make_graph(graph):
    
    ma_graph=nx.Graph()
    
    for patien in graph.columns:
        graph_0 = graph[graph[patien]!=0][patien]
        suma = np.sum(graph_0)
        s16S_norm=graph_0/suma
        
        g_m = nx.complete_graph(len(graph_0))
        
        labels={}
        for i, tax_i in zip(range(len(graph_0)),graph_0.index):
            labels[i]=tax_i
            g_m=nx.relabel_nodes(g_m,labels)
    
        for u,v,d in g_m.edges(data=True):
            w_uv = s16S_norm[u]*s16S_norm[v]
        
            if (u,v) not in ma_graph.edges():            
                ma_graph.add_edge(u,v,weight=w_uv)
            else:
                weight_h_old=ma_graph[u][v]['weight']                
                ma_graph[u][v].update(weight = weight_h_old + w_uv )
                
        #nx.write_gexf(ma_graph,'g_h_saliva.gexf')
    
    return ma_graph


# We write the network to grpah it in Gephi.
# df_h_hlung, df_h_saliva, df_s_saliva, df_s_hlung, df_s_alung, df_s_faeces
ma_graph = make_graph_median(df_s_faeces)
nx.write_gexf(ma_graph,'g_s_faeces_median.gexf')






