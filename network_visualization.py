import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import mysql.connector as connector
import matplotlib.pyplot as plt
import argparse
import os


my_hst = "localhost"
my_usr = "junting3"
my_pass = "FuEUylqI7aG2IAtf"
my_SGD_dbName = "Inparanoid_db"

def get_sc_dict(SC,connected):
    network_dict = {}
    for i in range(len(SC)):
	sc = SC[i]
	neighbor = connected[i]
	if sc not in network_dict.keys():
	    network_dict[sc] = []
	elif neighbor not in network_dict[sc]:
	    network_dict[sc].append(neighbor)
    return network_dict  


def get_adj_dict(network_dict):
    database_connection = connector.connect(host=my_hst,user=my_usr,passwd=my_pass)
    cursor = database_connection.cursor()
    cursor.execute("USE " + "SGD_db")
    adj_dict = {}
    for key in network_dict.keys():
	neighbors = network_dict[key]
	print("Length of neighbors ", len(neighbors))
	adj = np.zeros((len(neighbors),len(neighbors)))
	num_dict = {k: v for v, k in enumerate(neighbors)}
	for i in range(len(neighbors)):
	    curr_sc = neighbors[i]
	    cursor.execute("SELECT protein2 FROM string_network where protein1 LIKE '%"+ curr_sc +"%'")
	    result= cursor.fetchall()
	    for t in result:
		g = str(t[0].split(".")[1])
		if g in neighbors:
		    adj[i][num_dict[g]] = 1
		    adj[num_dict[g]][i] = 1
	adj_dict[key] = adj
    return adj_dict


def plot_adj(adj_dict,network_dict,map_dict,show_full,path):
    map_path = path
    if show_full:
        os.mkdir(path+"mapping/")
        map_path = path+"mapping/"
    for key in adj_dict.keys():
	if show_full:
            f = open(map_path+key+".txt","w+")
            f.write("Target"+"\t" + "SC"+"\t"+"Connected_SC"+'\n')
	adj = adj_dict[key]
	neighbors = network_dict[key]
	num_dict = {k: v for v, k in enumerate(neighbors)}
	num_dict = {v: k for k, v in num_dict.iteritems()}
	rowsum = np.sum(adj_dict[key],axis=1)
	ind = range(len(rowsum))
	color_map = []
	label_dict = {}
	if len(rowsum) > 10:
	    ind = np.argpartition(rowsum, len(rowsum) - 10)[-10:]
	G = nx.Graph()
	G.add_node(key)
	label_dict[key] = key
	nb = []
	for i in ind:
	    #G.add_node(num_dict[i]+"\n"+ str(1))
	    nb.append(num_dict[i])
	    
	    if len(map_dict[key][num_dict[i]])>1:
	        label_dict [num_dict[i]] = num_dict[i]+"\n"+ str(map_dict[key][num_dict[i]][0])+ "\n"+ str(map_dict[key][num_dict[i]][1])
	    else: 
	        label_dict [num_dict[i]] = num_dict[i]+"\n"+ str(map_dict[key][num_dict[i]][0])
	    if show_full:
		for j in map_dict[key][num_dict[i]]:
		    f.write(j+"\t" +key+"\t" +num_dict[i]+'\n')
	for i in ind:
	    G.add_edge(key,num_dict[i])
	    for j in ind:
		if adj[i][j] ==1:
		    G.add_edge(num_dict[i],num_dict[j]) 
        pos = nx.spring_layout(G,k=0.95)
	nx.draw_networkx_edges(G,pos,width=0.3)
	nx.draw_networkx_nodes(G,pos,nodelist=[key],node_color='b', node_size=800)
	nx.draw_networkx_nodes(G,pos,nodelist=nb,node_color='r', node_size=1000)
	nx.draw_networkx_labels(G,pos, labels = label_dict, font_weight = "light",font_size=6)
	plt.axis("off")
	plt.savefig(path+str(key)+".pdf")
	plt.clf()


def read_file(query):
    f = open("fullresult_"+query+".txt",'r') 
    lines = f.readlines()
    SC = []
    connected = []
    map_dict = {}
    for line in lines[1:]:
        line = line.split('\t')
	curr_sc = line[3]
	curr_connected = line[5]
        SC.append(curr_sc)
        connected.append(curr_connected)
	if curr_sc not in map_dict.keys():
	    map_dict[curr_sc] = {}
	    map_dict[curr_sc][curr_connected] = []
	    map_dict[curr_sc][curr_connected].append(line[2])
	else:
	    if curr_connected not in map_dict[curr_sc].keys():
		map_dict[curr_sc][curr_connected] = [line[2]]
	    else:
	        map_dict[curr_sc][curr_connected].append(line[2])
    return SC,connected,map_dict


if __name__ == "__main__":
 #   network_dict = get_sc_dict(SC,connected)
  #  adjancy_dict = get_adj_dict(network_dict)
   # np.save('/home/a-m/junting3/network_dict.npy', network_dict)
    #np.save('/home/a-m/junting3/adj_dict.npy', adjancy_dict) 
    parser = argparse.ArgumentParser()
    parser.add_argument('-show_full', dest='show_full', default=False, action='store_true')
    parser.add_argument("--query", type=str, default='205437')
    args = parser.parse_args()
    SC,connected,map_dict = read_file(args.query)   
    adjancy_dict = np.load("/home/a-m/junting3/adj_dict.npy").item()
    network_dict = np.load("/home/a-m/junting3/network_dict.npy").item()
    path = os.getcwd()
    path = path + "/"+args.query+"/"
    os.mkdir(path)
    plot_adj(adjancy_dict,network_dict,map_dict,args.show_full,path)
