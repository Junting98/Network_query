import mysql.connector as connector
import csv
import sys
import os
dict = {"0":"neighborhood", "1":"neighborhood_transferred", "2":"fusion", "3":"cooccurance", "4":"homology",
        "5":"coexpression","6":"coexpression_transferred","7":"experiments","8":"experiments_transferred","9":"database_origin",
        "10":"databse_transferred","11":"textmining","12":"textmining_transferred"}
table_dict = {0:"LP_SC",1:"RT_SC",2:"LP_genes", 3:"RT_genes"}

def query_network(hst,usr,password,SGD_dbName,query,table):
#### Get the SC orthologs
    
    database_connection = connector.connect(host=hst,user=usr,passwd=password)
    cursor = database_connection.cursor()
    cursor.execute("USE " + SGD_dbName)
    cursor.execute("SELECT SC_genes, Score_in FROM "+table_dict[table]+" where "+table_dict[table+2]+" LIKE '%"+query+"%'")
    result = cursor.fetchall() 
    SC_Orthologs = []
    SC_dict = {}
    for t in result:
	SC = str(t[0])
        SC_Orthologs.append(SC)
	SC_dict[SC] = t[1]
 
    return SC_Orthologs, SC_dict

def get_connected(hst,usr,password,SGD_dbName,SC_Orthologs,edges,threshold):
#####Seach for other connected SC
    database_connection = connector.connect(host=hst,user=usr,passwd=password)
    cursor = database_connection.cursor()
    cursor.execute("USE " + "SGD_db")
    edge_type = []
    #print(edges)
    if threshold is "":
	threshold = "0"
    if edges == [""]:
	edge_type.append("textmining")
    else:
        for edge in edges:
            edge_type.append(dict[edge])
    #print(edge_type)
    network = {}
    network_score = {}
    for edge in edge_type:
	#print "Current SC is " + SC
	network[edge] = {}
	network_score[edge] = {}
        for SC in SC_Orthologs:
	    cursor.execute("SELECT protein2, "+edge+" FROM string_network where protein1 LIKE '%"+SC+"%' AND "+edge+" > "+threshold)
	    result = cursor.fetchall()
            for t in result:
		g = str(t[0].split(".")[1])
		score = t[1]
		if SC in network[edge].keys():
		    network[edge][SC].append(g)
		    network_score[edge][SC].append(score)
		else:
		    network[edge][SC] = [g]
		    network_score[edge][SC] = [score]
    #print(network)
    #print "Number of connected Components" + str(len(network))
    #cursor.close()
    return network,edge_type,network_score


def get_query_orthologs(hst,usr,password,SGD_dbName,SC_Orthologs, connected_SC,table,edge_type ):
    database_connection = connector.connect(host=hst,user=usr,passwd=password)
    cursor = database_connection.cursor()
    cursor.execute("USE " + SGD_dbName)
    query_orthologs = {}
    query_score = {}
    for edge in edge_type:
	query_orthologs[edge] ={}
	query_score[edge] = {}
        for SC in connected_SC[edge].keys():
	    query_orthologs[edge][SC] = {}
	    query_score[edge][SC] = {}
	    for gene in connected_SC[edge][SC]:
		#print(gene)
	        cursor.execute("SELECT "+table_dict[table+2]+", Score_in FROM "+table_dict[table]+" where SC_genes LIKE '%"+gene+"%'")
	        result = cursor.fetchall()
	        if result != []:
	            for t in result:
			score = t[1]
			g = str(t[0].split("|")[2])
	        	if gene in query_orthologs[edge][SC] != {}:			    
			    query_orthologs[edge][SC][gene].append(g)
			    query_score[edge][SC][gene].append(score)
		        else:
			    query_orthologs[edge][SC][gene] = [g]
   			    query_score[edge][SC][gene] = [score] 
	        else:
		    query_orthologs[edge][SC][gene] = []
    #print "Number of orthologs in query species" + str(len(query_orthologs))
    #print(query_orthologs)
    return query_orthologs,query_score


def write_to_file(query_orthologs, query,SC_dict,SC_Orthologs, query_score,network,network_score):
    #print(query_orthologs)
    with open("fullresult_"+query+".txt",'w+') as f:
        f.write("Edge type" + "\t"+ "Query"+"\t" +"Target Species" + "\t"+ "SC" + "\t"+ "Scores for target Species to SC" + "\t"+ "Connected SC" + "\t"+ "Connection score" + "\t"+"Scores for SC to target Species" +"\n")
	for edge in network_score.keys():
	    for i in range(len(network[edge])):
		SC = SC_Orthologs[i]
		SC_score = SC_dict[SC]
		connectSC = network[edge][SC]
		connectScores = network_score[edge][SC]
		for j in range(len(connectSC)):
		    genes = connectSC[j]
		    genes_score = connectScores[j]
		    #print(SC)
		    #print(genes)
		    for k in range(len(query_orthologs[edge][SC][genes])):
			target = query_orthologs[edge][SC][genes][k]
			target_score = query_score[edge][SC][genes][k]
			f.write(edge + "\t"+ query + "\t"+target + "\t"+ SC + "\t"+ str(SC_score) + "\t"+ genes + "\t"+ str(genes_score) + "\t"+str(target_score) +"\n")
    f.close()


if __name__ == "__main__":
    my_hst = "localhost"
    my_usr = "junting3"
    my_pass = "FuEUylqI7aG2IAtf"
    my_SGD_dbName = "Inparanoid_db"
    create_db = False
    print("Which kind of edges to query:")
    print("Please enter the number:")
    print("0: neighborhood")
    print("1: neighborhood_transferred")
    print("2: fusion")
    print("3: cooccurance")
    print("4: homology")
    print("5: coexpression")
    print("6: coexpression_transferred")
    print("7: experiments")
    print("8: experiments_transferred")
    print("9: database_origin")
    print("10: databse_transferred")
    print("11: textmining")
    print("12: textmining_transferred")
    
    query = sys.argv[1]
    if create_db:
        create_database(my_hst, my_usr, my_pass, my_SGD_dbName)
    edges = raw_input("Please enter the number seperated by space. (default is 11: text mining)").strip("\n")
    threshold = raw_input("Please enter the threshold for the scores, 0 - 1000. (default is 0)" )
    table= raw_input("Which table to query: 0: LP_SC, 1:RT_SC. (default is 0)")
    if table is "":
	table = 0
    else: table = int(table)
    if " " in edges:
	edges = edges.split(" ")
    else:
	edges = [edges]

    database_connection = connector.connect(host=my_hst,user=my_usr,passwd=my_pass)
    cursor = database_connection.cursor()
    cursor.execute("USE " + my_SGD_dbName)
    cursor.execute("SELECT SC_genes, Score_in FROM "+table_dict[table]+" where "+table_dict[table+2]+" LIKE '%"+query+"%'")
    result = cursor.fetchall()
    if result == []:
	print("Currently no such gene exists in the data base.")
 	choice = raw_input("Do you want to run OrthoTator to get the result? (y/n)")
	print(choice)
	if choice.lower() == "n" or "no":
	    exit()
	command = raw_input("Please provide a command for running OrthoTator.")
	os.system(command)
    SC_Orthologs,SC_dict = query_network(my_hst, my_usr, my_pass, my_SGD_dbName,query,table)
    Connected_SC,edge_type,network_score =  get_connected(my_hst, my_usr, my_pass, my_SGD_dbName,SC_Orthologs,edges,threshold)
    query_orthologs,query_score = get_query_orthologs(my_hst, my_usr, my_pass, my_SGD_dbName,SC_Orthologs,Connected_SC,table,edge_type)
    write_to_file(query_orthologs, query,SC_dict,SC_Orthologs,query_score,Connected_SC,network_score)
   
