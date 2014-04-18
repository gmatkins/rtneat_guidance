package rtNEAT;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.Vector;

//#ifndef _GENOME_H_
//#define _GENOME_H_
//
//#include <vector>
//#include "gene.h"
//#include "innovation.h"
//
//namespace NEAT {
//
//	enum mutator {
//		GAUSSIAN = 0,
//		COLDGAUSSIAN = 1
//	};
//
//	//----------------------------------------------------------------------- 
//	//A Genome is the primary source of genotype information used to create   
//	//a phenotype.  It contains 3 major constituents:                         
//	//  1) A list of Traits                                                 
//	//  2) A list of NNodes pointing to a Trait from (1)                      
//	//  3) A list of Genes with Links that point to Traits from (1)           
//	//(1) Reserved parameter space for future use
//	//(2) NNode specifications                                                
//	//(3) Is the primary source of innovation in the evolutionary Genome.     
//	//    Each Gene in (3) has a marker telling when it arose historically.   
//	//    Thus, these Genes can be used to speciate the population, and the   
//	//    list of Genes provide an evolutionary history of innovation and     
//	//    link-building.
//
//	class Genome {
class Genome{
	
	public enum mutator {GAUSSIAN, COLDGAUSSIAN};
//
//	public:
//		int genome_id;
	public int genome_id;
//
//		std::vector<Trait*> traits; //parameter conglomerations
	public Vector<Trait> traits;
//		std::vector<NNode*> nodes; //List of NNodes for the Network
	public Vector<Nnode> nodes;
//		std::vector<Gene*> genes; //List of innovation-tracking genes
	public Vector<Gene> genes;
//
//		Network *phenotype; //Allows Genome to be matched with its Network
	public Network phenotype;
//
//		int get_last_node_id(); //Return id of final NNode in Genome
	public int get_last_node_id() {
		return nodes.lastElement().node_id;
	}
//		double get_last_gene_innovnum(); //Return last innovation number in Genome
	public double get_last_gene_innovnum() {
		//return ((*(genes.end() - 1))->innovation_num)+1;
		return genes.lastElement().innovation_num;
	}
//
//		void print_genome(); //Displays Genome on screen
//
//		//Constructor which takes full genome specs and puts them into the new one
//		Genome(int id, std::vector<Trait*> t, std::vector<NNode*> n, std::vector<Gene*> g);
//	public Genome(int id, Vector<Trait> t, Vector<Nnode> n, Vector<Gene> g) {
//		genome_id=id;
//		traits= new Vector<Trait>(t);
//		nodes= new Vector<Nnode>(n); 
//		genes= new Vector<Gene>(g);
//	}
//
//		//Constructor which takes in links (not genes) and creates a Genome
//		Genome(int id, std::vector<Trait*> t, std::vector<NNode*> n, std::vector<Link*> links);
	public Genome(int id, Vector<Trait> t, Vector<Nnode> n, Vector<Gene> g, Vector<Link> links) {
		//std::vector<Link*>::iterator curlink;
		Gene tempgene;
		if (t != null) traits= new Vector<Trait>(t);
		if (n != null) nodes= new Vector<Nnode>(n); 
		if (g != null) genes= new Vector<Gene>(g);

		genome_id=id;

		//We go through the links and turn them into original genes
		//for(curlink=links.begin();curlink!=links.end();++curlink) {
		for (Link curlink : links){
			//Create genes one at a time
			tempgene=new Gene((curlink).linktrait, (curlink).weight,(curlink).in_node,(curlink).out_node,(curlink).is_recurrent,1.0,0.0);
			genes.add(tempgene);
		}

	}
//
//		// Copy constructor
//		Genome(const Genome& genome);
	
	public Genome(Genome genome)
	{
		genome_id = genome.genome_id;

//		Vector<Trait>::const_iterator curtrait;
//		std::vector<NNode*>::const_iterator curnode;
//		std::vector<Gene*>::const_iterator curgene;

		//for(curtrait=genome.traits.begin(); curtrait!=genome.traits.end(); ++curtrait) {
		for (Trait curtrait : genome.traits){
			traits.add(new Trait(curtrait));
		}

		Trait assoc_trait;
		//Duplicate NNodes
		//for(curnode=genome.nodes.begin();curnode!=genome.nodes.end();++curnode) {
		for(Nnode curnode : genome.nodes){
			//First, find the trait that this node points to
			if (((curnode).nodetrait)==null) assoc_trait=null;
			else {
				//curtrait=traits.firstElement();
//				while(((*curtrait)->trait_id)!=(((*curnode)->nodetrait)->trait_id))
//					++curtrait;
				for(Trait curtrait : genome.traits){
					if (curtrait.trait_id != curnode.nodetrait.trait_id) continue;
					assoc_trait = new Trait(curtrait);
					break;
				}
				//assoc_trait=(*curtrait);
			}

			Nnode newnode=new Nnode(curnode,assoc_trait);

			(curnode).dup=newnode;  //Remember this node's old copy
			//    (*curnode)->activation_count=55;
			nodes.add(newnode);    
		}

		Nnode inode; //For forming a gene 
		Nnode onode; //For forming a gene
		Trait traitptr;

		//Duplicate Genes
		//for(curgene=genome.genes.begin(); curgene!=genome.genes.end(); ++curgene) {
		for (Gene curgene : genome.genes){
			//First find the nodes connected by the gene's link

			inode=(((curgene).lnk).in_node).dup;
			onode=(((curgene).lnk).out_node).dup;

			//Get a pointer to the trait expressed by this gene
			traitptr=((curgene).lnk).linktrait;
			if (traitptr==null) assoc_trait=null;
			else {
//				curtrait=traits.begin();
//				while(((*curtrait)->trait_id)!=(traitptr->trait_id))
//					++curtrait;
//				assoc_trait=(*curtrait);
				for (Trait curtrait : traits){
					if (curtrait.trait_id != traitptr.trait_id) continue;
					assoc_trait = new Trait(curtrait);
					break;
				}
			}

			Gene newgene=new Gene(curgene,assoc_trait,inode,onode);
			genes.add(newgene);

		}
	}

//
//		//Special constructor which spawns off an input file
//		//This constructor assumes that some routine has already read in GENOMESTART
//        Genome(int id, std::ifstream &iFile);
	public Genome(int id, BufferedReader iFile) {

		//char curword[128];  //max word size of 128 characters
		String curword;
//		char curline[1024]; //max line size of 1024 characters
		String curline;
//		char delimiters[] = " \n";
		String delimeters = new String(" \n");

		//int done=0;
		boolean done = false;

		//int pause;

		genome_id=id;

		//iFile.getline(curline, sizeof(curline));
		curline = iFile.readLine();
		//int wordcount = NEAT::getUnitCount(curline, delimiters);
		int wordcount = Neat.getUnitCount(curline, delimeters);
		int curwordnum = 0;

		//Loop until file is finished, parsing each line
		while (!done) {

	        //std::cout << curline << std::endl;

			if (curwordnum > wordcount || wordcount == 0) {
				//iFile.getline(curline, sizeof(curline));
				curline = iFile.readLine();
				//wordcount = NEAT::getUnitCount(curline, delimiters);
				wordcount = Neat.getUnitCount(curline, delimeters);
				curwordnum = 0;
			}
	        
	        //std::stringstream ss(curline);
			//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
	       // ss >> curword;
			curword = curline.split(" ", 2)[0];
			curline = new String(curline.split(" ", 2)[1]);

			//printf(curword);
			//printf(" test\n");
			//Check for end of Genome
			//if (strcmp(curword,"genomeend")==0) {
			if (curword == "genomeend"){
				//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
	            //ss >> curword;
				curword = curline.split(" ", 2)[0];
				curline = new String(curline.split(" ", 2)[1]);
				int idcheck = Integer.parseInt(curword);
				//iFile>>idcheck;
				if (idcheck!=genome_id) System.out.println("ERROR: id mismatch in genome");
				done=true;
			}

			//Ignore genomestart if it hasn't been gobbled yet
			else if ((curword == "genomestart")) {
				++curwordnum;
				//cout<<"genomestart"<<endl;
			}

			//Ignore comments surrounded by - they get printed to screen
			else if ((curword == "/*")) {
				//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
	            //ss >> curword;
	            curword = curline.split(" ", 2)[0];
				curline = new String(curline.split(" ", 2)[1]);
				while ((curword != "*/")) {
					//cout<<curword<<" ";
					//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
	                //ss >> curword;
	                curword = curline.split(" ", 2)[0];
					curline = new String(curline.split(" ", 2)[1]);
				}
				//cout<<endl;
			}

			//Read in a trait
			else if ((curword == "trait")) {
				Trait newtrait;

				//char argline[1024];
				String argline;
				//strcpy(argline, NEAT::getUnits(curline, curwordnum, wordcount, delimiters));

				curwordnum = wordcount + 1;

	            //ss.getline(argline, 1024);
				argline = curline.split("\n", 2)[0];
				curline = curline.split("\n", 2)[1];
				//Allocate the new trait
				newtrait=new Trait(argline);

				//Add trait to vector of traits
				traits.add(newtrait);
			}

			//Read in a node
			else if ((curword == "node")) {
				Nnode newnode;

				//char argline[1024];
				String argline;
				//strcpy(argline, NEAT::getUnits(curline, curwordnum, wordcount, delimiters));
				curwordnum = wordcount + 1;
	            
	            //ss.getline(argline, 1024);
				argline = curline.split("\n", 2)[0];
				curline = curline.split("\n", 2)[1];
				//Allocate the new node
				newnode=new Nnode(argline,traits);

				//Add the node to the list of nodes
				nodes.add(newnode);
			}

			//Read in a Gene
			else if ((curword == "gene")) {
				Gene newgene;

				//char argline[1024];
				String argline;
				//strcpy(argline, NEAT::getUnits(curline, curwordnum, wordcount, delimiters));
				curwordnum = wordcount + 1;

	            //ss.getline(argline, 1024);
				argline = curline.split("\n", 2)[0];
				curline = curline.split("\n", 2)[1];
	            //std::cout << "New gene: " << ss.str() << std::endl;
				//Allocate the new Gene
	            newgene=new Gene(argline,traits,nodes);

				//Add the gene to the genome
				genes.add(newgene);

	            //std::cout<<"Added gene " << newgene << std::endl;
			}

		}

	}
//
//		// This special constructor creates a Genome
//		// with i inputs, o outputs, n out of nmax hidden units, and random
//		// connectivity.  If r is true then recurrent connections will
//		// be included. 
//		// The last input is a bias
//		// Linkprob is the probability of a link  
//		Genome(int new_id,int i, int o, int n,int nmax, bool r, double linkprob);
	public Genome(int new_id,int i, int o, int n,int nmax, bool r, double linkprob) {
		int totalnodes;
		boolean[] cm; //The connection matrix which will be randomized
		boolean[] cmp; //Connection matrix pointer
		int matrixdim;
		int count;

		int ncount; //Node and connection counters
		int ccount;

		int row;  //For navigating the matrix
		int col;

		double new_weight;

		int maxnode; //No nodes above this number for this genome

		int first_output; //Number of first output node

		totalnodes=i+o+nmax;
		matrixdim=totalnodes*totalnodes;
		cm=new boolean[matrixdim];  //Dimension the connection matrix
		maxnode=i+n;

		first_output=totalnodes-o+1;

		//For creating the new genes
		Nnode newnode;
		Gene newgene;
		Trait newtrait;
		NNode in_node;
		NNode out_node;

		//Retrieves the nodes pointed to by connection genes
		//Vector<Nnode>::iterator node_iter;

		//Assign the id
		genome_id=new_id;

		//cout<<"Assigned id "<<genome_id<<endl;

		//Step through the connection matrix, randomly assigning bits
		cmp=cm;
		for(count=0;count<matrixdim;count++) {
			if (Neat.randfloat()<linkprob)
				cmp[count]=true;
			else cmp[count]=false;
			//cmp++;
		}

		//Create a dummy trait (this is for future expansion of the system)
		newtrait=new Trait(1,0,0,0,0,0,0,0,0,0);
		traits.add(newtrait);

		//Build the input nodes
		for(ncount=1;ncount<=i;ncount++) {
			if (ncount<i)
				newnode=new Nnode(Nnode.nodetype.SENSOR,ncount,Nnode.nodeplace.INPUT);
			else newnode=new Nnode(Nnode.nodetype.SENSOR,ncount,Nnode.nodeplace.BIAS);

			newnode.nodetrait=newtrait;

			//Add the node to the list of nodes
			nodes.add(newnode);
		}

		//Build the hidden nodes
		for(ncount=i+1;ncount<=i+n;ncount++) {
			newnode=new Nnode(Nnode.nodetype.NEURON,ncount,Nnode.nodeplace.HIDDEN);
			newnode.nodetrait=newtrait;
			//Add the node to the list of nodes
			nodes.add(newnode);
		}

		//Build the output nodes
		for(ncount=first_output;ncount<=totalnodes;ncount++) {
			newnode=new Nnode(Nnode.nodetype.NEURON,ncount,Nnode.nodeplace.OUTPUT);
			newnode.nodetrait=newtrait;
			//Add the node to the list of nodes
			nodes.add(newnode);
		}

		//cout<<"Built nodes"<<endl;

		//Connect the nodes 
		ccount=1;  //Start the connection counter

		//Step through the connection matrix, creating connection genes
		cmp=cm;
		count=0;
		for(col=1;col<=totalnodes;col++)
			for(row=1;row<=totalnodes;row++) {
				//Only try to create a link if it is in the matrix
				//and not leading into a sensor

				if ((cmp[row]==true)&&(col>i)&&
					((col<=maxnode)||(col>=first_output))&&
					((row<=maxnode)||(row>=first_output))) {
						//If it isn't recurrent, create the connection no matter what
						if (col>row) {

							//Retrieve the in_node
//							node_iter=nodes.begin();
//							while((*node_iter)->node_id!=row)
//								node_iter++;
//
//							in_node=(*node_iter);
							for (Nnode node_iter : nodes){
								if(node_iter.node_id != row) continue;
								in_node = new Nnode(node_iter);
								break;
							}

							//Retrieve the out_node
//							node_iter=nodes.begin();
//							while((*node_iter)->node_id!=col)
//								node_iter++;
//
//							out_node=(*node_iter);
							for (Nnode node_iter : nodes){
								if(node_iter.node_id != col) continue;
								out_node = new Nnode(node_iter);
								break;
							}

							//Create the gene
							new_weight=Neat.randposneg() * Neat.randfloat();
							newgene=new Gene(newtrait,new_weight, in_node, out_node,false,count,new_weight);

							//Add the gene to the genome
							genes.add(newgene);
						}
						else if (r) {
							//Create a recurrent connection

							//Retrieve the in_node
//							node_iter=nodes.begin();
//							while((*node_iter)->node_id!=row)
//								node_iter++;
//
//							in_node=(*node_iter);
							for (Nnode node_iter : nodes){
								if(node_iter.node_id != row) continue;
								in_node = new Nnode(node_iter);
								break;
							}

							//Retrieve the out_node
//							node_iter=nodes.begin();
//							while((*node_iter)->node_id!=col)
//								node_iter++;
//
//							out_node=(*node_iter);
							for (Nnode node_iter : nodes){
								if(node_iter.node_id != col) continue;
								out_node = new Nnode(node_iter);
								break;
							}

							//Create the gene
							new_weight=Neat.randposneg() * Neat.randfloat();
							newgene=new Gene(newtrait,new_weight, in_node, out_node,true,count,new_weight);

							//Add the gene to the genome
							genes.add(newgene);

						}

					}

					count++; //increment gene counter	    
					//cmp++;
			}

			//delete [] cm;

	}
//
//		//Special constructor that creates a Genome of 3 possible types:
//		//0 - Fully linked, no hidden nodes
//		//1 - Fully linked, one hidden node splitting each link
//		//2 - Fully connected with a hidden layer, recurrent 
//		//num_hidden is only used in type 2
//		Genome(int num_in,int num_out,int num_hidden,int type);
	public Genome(int num_in,int num_out,int num_hidden,int type) {

		//Temporary lists of nodes
		Vector<Nnode> inputs;
		Vector<Nnode> outputs;
		Vector<Nnode> hidden;
		Nnode bias; //Remember the bias

//		std::vector<NNode*>::iterator curnode1; //Node iterator1
//		std::vector<NNode*>::iterator curnode2; //Node iterator2
//		std::vector<NNode*>::iterator curnode3; //Node iterator3

		//For creating the new genes
		Nnode newnode;
		Gene newgene;
		Trait newtrait;

		int count;
		int ncount;


		//Assign the id 0
		genome_id=0;

		//Create a dummy trait (this is for future expansion of the system)
		newtrait=new Trait(1,0,0,0,0,0,0,0,0,0);
		traits.add(newtrait);

		//Adjust hidden number
		if (type==0) 
			num_hidden=0;
		else if (type==1)
			num_hidden=num_in*num_out;

		//Create the inputs and outputs

		//Build the input nodes
		for(ncount=1;ncount<=num_in;ncount++) {
			if (ncount<num_in)
				newnode=new Nnode(Nnode.nodetype.SENSOR,ncount,Nnode.nodeplace.INPUT);
			else { 
				newnode=new Nnode(Nnode.nodetype.SENSOR,ncount,Nnode.nodeplace.BIAS);
				bias=newnode;
			}

			//newnode->nodetrait=newtrait;

			//Add the node to the list of nodes
			nodes.add(newnode);
			inputs.add(newnode);
		}

		//Build the hidden nodes
		for(ncount=num_in+1;ncount<=num_in+num_hidden;ncount++) {
			newnode=new Nnode(Nnode.nodetype.NEURON,ncount,Nnode.nodeplace.HIDDEN);
			//newnode->nodetrait=newtrait;
			//Add the node to the list of nodes
			nodes.add(newnode);
			hidden.add(newnode);
		}

		//Build the output nodes
		for(ncount=num_in+num_hidden+1;ncount<=num_in+num_hidden+num_out;ncount++) {
			newnode=new Nnode(Nnode.nodetype.NEURON,ncount,Nnode.nodeplace.OUTPUT);
			//newnode->nodetrait=newtrait;
			//Add the node to the list of nodes
			nodes.add(newnode);
			outputs.add(newnode);
		}

		//Create the links depending on the type
		if (type==0) {
			//Just connect inputs straight to outputs

			count=1;

			//Loop over the outputs
			//for(curnode1=outputs.begin();curnode1!=outputs.end();++curnode1) {
			for (Nnode curnode1 : outputs){
				//Loop over the inputs
				//for(curnode2=inputs.begin();curnode2!=inputs.end();++curnode2) {
				for (Nnode curnode2 : inputs){
					//Connect each input to each output
					newgene=new Gene(newtrait,0, (curnode2), (curnode1),false,count,0);

					//Add the gene to the genome
					genes.add(newgene);	 

					count++;

				}

			}

		} //end type 0
		//A split link from each input to each output
		else if (type==1) {
			count=1; //Start the gene number counter

			//curnode3=hidden.begin(); //One hidden for ever input-output pair
			int curnode3 = 0;
			//Loop over the outputs
			//for(curnode1=outputs.begin();curnode1!=outputs.end();++curnode1) {
			for (Nnode curnode1 : outputs){
				//Loop over the inputs
				//for(curnode2=inputs.begin();curnode2!=inputs.end();++curnode2) {
				for (Nnode curnode2 : inputs){

					//Connect Input to hidden
					newgene=new Gene(newtrait,0, (curnode2), (hidden.get(curnode3)),false,count,0);
					//Add the gene to the genome
					genes.add(newgene);

					count++; //Next gene

					//Connect hidden to output
					newgene=new Gene(newtrait,0, (hidden.get(curnode3)), (curnode1),false,count,0);
					//Add the gene to the genome
					genes.add(newgene);

					++curnode3; //Next hidden node
					count++; //Next gene

				}
			}

		}//end type 1
		//Fully connected 
		else if (type==2) {
			count=1; //Start gene counter at 1


			//Connect all inputs to all hidden nodes
			//for(curnode1=hidden.begin();curnode1!=hidden.end();++curnode1) {
			for (Nnode curnode1 : hidden){
				//Loop over the inputs
				//for(curnode2=inputs.begin();curnode2!=inputs.end();++curnode2) {
				for(Nnode curnode2 : inputs){
					//Connect each input to each hidden
					newgene=new Gene(newtrait,0, (curnode2), (curnode1),false,count,0);

					//Add the gene to the genome
					genes.add(newgene);	 

					count++;

				}
			}

			//Connect all hidden units to all outputs
			//for(curnode1=outputs.begin();curnode1!=outputs.end();++curnode1) {
			for (Nnode curnode1 ; outputs){
				//Loop over the inputs
				//for(curnode2=hidden.begin();curnode2!=hidden.end();++curnode2) {
				for (Nnode curnode2 : hidden){
					//Connect each input to each hidden
					newgene=new Gene(newtrait,0, (curnode2), (curnode1),false,count,0);

					//Add the gene to the genome
					genes.add(newgene);	 

					count++;

				}
			}

			//Connect the bias to all outputs
			//for(curnode1=outputs.begin();curnode1!=outputs.end();++curnode1) {
			for (Nnode curnode1 : outputs){
				newgene=new Gene(newtrait,0, bias, (curnode1),false,count,0);

				//Add the gene to the genome
				genes.add(newgene);	 

				count++;
			}

			//Recurrently connect the hidden nodes
			//for(curnode1=hidden.begin();curnode1!=hidden.end();++curnode1) {
			for (Nnode curnode1 : hidden){
				//Loop Over all Hidden
				//for(curnode2=hidden.begin();curnode2!=hidden.end();++curnode2) {
				for (Nnode curnode2 : hidden){
					//Connect each hidden to each hidden
					newgene=new Gene(newtrait,0, (curnode2), (curnode1),true,count,0);

					//Add the gene to the genome
					genes.add(newgene);	 

					count++;

				}

			}

		}//end type 2

	}
//
//		// Loads a new Genome from a file (doesn't require knowledge of Genome's id)
//		static Genome *new_Genome_load(char *filename);
	public static Genome new_Genome_load(String filename) {
		Genome newgenome;

		int id;

		//char curline[1024];
		char curword[] = new char[20];  //max word size of 20 characters
		//String curword;
		//char delimiters[] = " \n";
		//int curwordnum = 0;

		//std::ifstream iFile(filename);
		try{
			InputStream is = new FileInputStream("filename");
			BufferedReader iFile = new BufferedReader(new InputStreamReader(is));

		//Make sure it worked
		//if (!iFile) {
		//	cerr<<"Can't open "<<filename<<" for input"<<endl;
		//	return 0;
		//}

		//iFile>>curword;
		iFile.read(curword);	
		//iFile.getline(curline, sizeof(curline));
		//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));

		//Bypass initial comment
		if ((new String(curword) == "/*")) {
			//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
			//iFile>>curword;
			iFile.read(curword);
			while ((new String(curword) != "*/")) {
				System.out.println(curword);
				//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
				//iFile>>curword;
				iFile.read(curword);
			}

			//cout<<endl;
			//iFile>>curword;
			iFile.read(curword);
			//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
		}

		//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
		//id = atoi(curword);
		//iFile>>id;
		id = iFile.read();

		newgenome=new Genome(id,iFile);

		iFile.close();
		return newgenome;
		}catch(Exception e){
			System.out.println("trouble constructing new Genome from file");
			System.out.println(e.getMessage());
			return null;
		}	
	}
//
//		//Destructor kills off all lists (including the trait vector)
//		~Genome();
//
//		//Generate a network phenotype from this Genome with specified id
//		Network *genesis(int);
	public Network genesis(int id) {
		//std::vector<NNode*>::iterator curnode; 
		//std::vector<Gene*>::iterator curgene;
		Nnode newnode;
		Trait curtrait;
		Link curlink;
		Link newlink;

		double maxweight=0.0; //Compute the maximum weight for adaptation purposes
		double weight_mag; //Measures absolute value of weights

		//Inputs and outputs will be collected here for the network
		//All nodes are collected in an all_list- 
		//this will be used for later safe destruction of the net
		Vector<Nnode> inlist;
		Vector<Nnode> outlist;
		Vector<Nnode> all_list;

		//Gene translation variables
		Nnode inode;
		Nnode onode;

		//The new network
		Network newnet;

		//Create the nodes
		//for(curnode=nodes.begin();curnode!=nodes.end();++curnode) {
		for(Nnode curnode : nodes){
			newnode=new Nnode((curnode).type,(curnode).node_id);

			//Derive the node parameters from the trait pointed to
			curtrait=(curnode).nodetrait;
			newnode.derive_trait(curtrait);

			//Check for input or output designation of node
			if (((curnode).gen_node_label)== Nnode.nodeplace.INPUT) 
				inlist.add(newnode);
			if (((curnode).gen_node_label)== Nnode.nodeplace.BIAS) 
				inlist.add(newnode);
			if (((curnode).gen_node_label)== Nnode.nodeplace.OUTPUT)
				outlist.add(newnode);

			//Keep track of all nodes, not just input and output
			all_list.add(newnode);

			//Have the node specifier point to the node it generated
			(curnode).analogue=newnode;

		}

		//Create the links by iterating through the genes
		//for(curgene=genes.begin();curgene!=genes.end();++curgene) {
		for (Gene curgene : genes){
			//Only create the link if the gene is enabled
			if (((curgene).enable)==true) {
				curlink=(curgene).lnk;
				inode=(curlink.in_node).analogue;
				onode=(curlink.out_node).analogue;
				//NOTE: This line could be run through a recurrency check if desired
				// (no need to in the current implementation of NEAT)
				newlink=new Link(curlink.weight,inode,onode,curlink.is_recurrent);

				(onode.incoming).add(newlink);
				(inode.outgoing).add(newlink);

				//Derive link's parameters from its Trait pointer
				curtrait=(curlink.linktrait);

				newlink.derive_trait(curtrait);

				//Keep track of maximum weight
				if (newlink.weight > 0)
					weight_mag=newlink.weight;
				else weight_mag=-newlink.weight;
				if (weight_mag>maxweight)
					maxweight=weight_mag;
			}
		}

		//Create the new network
		newnet=new Network(inlist,outlist,all_list,id);

		//Attach genotype and phenotype together
		newnet.genotype=this;
		phenotype=newnet;

		newnet.maxweight=maxweight;

		return newnet;

	}
//
//		// Dump this genome to specified file
//		void print_to_file(std::ostream &outFile);
//		void print_to_file(std::ofstream &outFile);
	public void print_to_file(BufferedWriter outFile) {
		  //std::vector<Trait*>::iterator curtrait;
		  //std::vector<NNode*>::iterator curnode;
		  //std::vector<Gene*>::iterator curgene;

		  //outFile<<"genomestart "<<genome_id<<std::endl;
		outFile.write("genomestart " + genome_id + "\n");

		  //Output the traits
		  //for(curtrait=traits.begin();curtrait!=traits.end();++curtrait) {
		for(Trait curtrait : traits){
//		    (curtrait).trait_id=curtrait-traits.begin()+1;
//		    (*curtrait)->print_to_file(outFile);
			curtrait.print_to_file(outFile);
		  }

		  //Output the nodes
		  //for(curnode=nodes.begin();curnode!=nodes.end();++curnode) {
		for(Nnode curnode : nodes){
		    //(*curnode)->print_to_file(outFile);
			curnode.print_to_file(outFile);
		  }

		  //Output the genes
		  //for(curgene=genes.begin();curgene!=genes.end();++curgene) {
		for(Gene curgene : genes){
		    //(*curgene)->print_to_file(outFile);
			curgene.print_to_file(outFile);
		  }

		  //outFile<<"genomeend "<<genome_id<<std::endl;
		outFile.write("genomeend" + genome_id + "\n");

		}
//
//		// Wrapper for print_to_file above
//		void print_to_filename(char *filename);
	public void print_to_filename(String filename) {
		//std::ofstream oFile(filename);
		//oFile.open(filename, std::ostream::Write);
		try{
			OutputStream is = new FileOutputStream("filename");
			BufferedWriter oFile = new BufferedWriter(new OutputStreamWriter(is));
			print_to_file(oFile);
			oFile.close();
		}catch(Exception e){
			System.out.println("Trouble writing to output file in Genome");
			System.out.println(e.getMessage());
		}
	}
//
//		// Duplicate this Genome to create a new one with the specified id 
//		Genome *duplicate(int new_id);
//
//		// For debugging: A number of tests can be run on a genome to check its
//		// integrity
//		// Note: Some of these tests do not indicate a bug, but rather are meant
//		// to be used to detect specific system states
//		bool verify();
	public boolean verify() {
		//std::vector<NNode*>::iterator curnode;
		//std::vector<Gene*>::iterator curgene;
		//std::vector<Gene*>::iterator curgene2;
		Nnode inode;
		Nnode onode;

		boolean disab;

		int last_id;

		//int pause;

		//cout<<"Verifying Genome id: "<<this->genome_id<<endl;

		if (this==null) {
			//cout<<"ERROR GENOME EMPTY"<<endl;
			//cin>>pause;
		}

		//Check each gene's nodes
		//for(curgene=genes.begin();curgene!=genes.end();++curgene) {
		for (Gene curgene : genes){
			inode=((curgene).lnk).in_node;
			onode=((curgene).lnk).out_node;

			//Look for inode
//			curnode=nodes.begin();
//			while((curnode!=nodes.end())&&
//				((*curnode)!=inode))
//				++curnode;
//
//			if (curnode==nodes.end()) {
//				//cout<<"MISSING iNODE FROM GENE NOT IN NODES OF GENOME!!"<<endl;
//				//cin>>pause;
//				return false;
//			}
			for (Nnode curnode : nodes){
				if (curnode!=nodes.lastElement() && curnode != inode) continue;
				if (curnode == nodes.lastElement()){
					return false;
				}
			}

			//Look for onode
//			curnode=nodes.begin();
//			while((curnode!=nodes.end())&&
//				((*curnode)!=onode))
//				++curnode;
//
//			if (curnode==nodes.end()) {
//				//cout<<"MISSING oNODE FROM GENE NOT IN NODES OF GENOME!!"<<endl;
//				//cin>>pause;
//				return false;
//			}
			for (Nnode curnode : nodes){
				if (curnode!=nodes.lastElement() && curnode != onode) continue;
				if (curnode == nodes.lastElement()){
					return false;
				}
			}

		}

		//Check for NNodes being out of order
		last_id=0;
		//for(curnode=nodes.begin();curnode!=nodes.end();++curnode) {
		for(Nnode curnode : nodes){
			if ((curnode).node_id<last_id) {
				//cout<<"ALERT: NODES OUT OF ORDER in "<<this<<endl;
				//cin>>pause;
				return false;
			}

			last_id=(curnode).node_id;
		}


		//Make sure there are no duplicate genes
		//for(curgene=genes.begin();curgene!=genes.end();++curgene) {
		for (Gene curgene : genes){

			//for(curgene2=genes.begin();curgene2!=genes.end();++curgene2) {
			for (Gene curgene2 : genes){
				if (((curgene)!=(curgene2))&&
					((((curgene).lnk).is_recurrent)==(((curgene2).lnk).is_recurrent))&&
					((((((curgene2).lnk).in_node).node_id)==((((curgene).lnk).in_node).node_id))&&
					(((((curgene2).lnk).out_node).node_id)==((((curgene).lnk).out_node).node_id)))) {
						//cout<<"ALERT: DUPLICATE GENES: "<<(*curgene)<<(*curgene2)<<endl;
						//cout<<"INSIDE GENOME: "<<this<<endl;

						//cin>>pause;
					}


			}
		}

		//See if a gene is not disabled properly
		//Note this check does not necessary mean anything is wrong
		//
		//if (nodes.size()>=15) {
		//disab=false;
		////Go through genes and see if one is disabled
		//for(curgene=genes.begin();curgene!=genes.end();++curgene) {
		//if (((*curgene)->enable)==false) disab=true;
		//}

		//if (disab==false) {
		//cout<<"ALERT: NO DISABLED GENE IN GENOME: "<<this<<endl;
		////cin>>pause;
		//}

		//}
		//

		//Check for 2 disables in a row
		//Note:  Again, this is not necessarily a bad sign
		if (nodes.size()>=500) {
			disab=false;
			//for(curgene=genes.begin();curgene!=genes.end();++curgene) {
			for (Gene curgene : genes){
				if ((((curgene).enable)==false)&&(disab==true)) {
					//cout<<"ALERT: 2 DISABLES IN A ROW: "<<this<<endl;
				}
				if (((curgene).enable)==false) disab=true;
				else disab=false;
			}
		}

		//cout<<"GENOME OK!"<<endl;

		return true;
	}
//
//		// ******* MUTATORS *******
//
//		// Perturb params in one trait
//		void mutate_random_trait();
//
//		// Change random link's trait. Repeat times times
//		void mutate_link_trait(int times);
//
//		// Change random node's trait times times 
//		void mutate_node_trait(int times);
//
//		// Add Gaussian noise to linkweights either GAUSSIAN or COLDGAUSSIAN (from zero)
//		void mutate_link_weights(double power,double rate,mutator mut_type);
//
//		// toggle genes on or off 
//		void mutate_toggle_enable(int times);
//
//		// Find first disabled gene and enable it 
//		void mutate_gene_reenable();
//
//		// These last kinds of mutations return false if they fail
//		//   They can fail under certain conditions,  being unable
//		//   to find a suitable place to make the mutation.
//		//   Generally, if they fail, they can be called again if desired. 
//
//		// Mutate genome by adding a node respresentation 
//		bool mutate_add_node(std::vector<Innovation*> &innovs,int &curnode_id,double &curinnov);
//
//		// Mutate the genome by adding a new link between 2 random NNodes 
//		bool mutate_add_link(std::vector<Innovation*> &innovs,double &curinnov,int tries); 
//
//		void mutate_add_sensor(std::vector<Innovation*> &innovs, double &curinnov);
//
//		// ****** MATING METHODS ***** 
//
//		// This method mates this Genome with another Genome g.  
//		//   For every point in each Genome, where each Genome shares
//		//   the innovation number, the Gene is chosen randomly from 
//		//   either parent.  If one parent has an innovation absent in 
//		//   the other, the baby will inherit the innovation 
//		//   Interspecies mating leads to all genes being inherited.
//		//   Otherwise, excess genes come from most fit parent.
//		Genome *mate_multipoint(Genome *g,int genomeid,double fitness1, double fitness2, bool interspec_flag);
//
//		//This method mates like multipoint but instead of selecting one
//		//   or the other when the innovation numbers match, it averages their
//		//   weights 
//		Genome *mate_multipoint_avg(Genome *g,int genomeid,double fitness1,double fitness2, bool interspec_flag);
//
//		// This method is similar to a standard single point CROSSOVER
//		//   operator.  Traits are averaged as in the previous 2 mating
//		//   methods.  A point is chosen in the smaller Genome for crossing
//		//   with the bigger one.  
//		Genome *mate_singlepoint(Genome *g,int genomeid);
//
//
//		// ******** COMPATIBILITY CHECKING METHODS ********
//
//		// This function gives a measure of compatibility between
//		//   two Genomes by computing a linear combination of 3
//		//   characterizing variables of their compatibilty.
//		//   The 3 variables represent PERCENT DISJOINT GENES, 
//		//   PERCENT EXCESS GENES, MUTATIONAL DIFFERENCE WITHIN
//		//   MATCHING GENES.  So the formula for compatibility 
//		//   is:  disjoint_coeff*pdg+excess_coeff*peg+mutdiff_coeff*mdmg.
//		//   The 3 coefficients are global system parameters 
//		double compatibility(Genome *g);
//
//		double trait_compare(Trait *t1,Trait *t2);
//
//		// Return number of non-disabled genes 
//		int extrons();
//
//		// Randomize the trait pointers of all the node and connection genes 
//		void randomize_traits();
//
//	protected:
//		//Inserts a NNode into a given ordered list of NNodes in order
//		void node_insert(std::vector<NNode*> &nlist, NNode *n);
//
//		//Adds a new gene that has been created through a mutation in the
//		//*correct order* into the list of genes in the genome
//		void add_gene(std::vector<Gene*> &glist,Gene *g);
//
//	};
	
} // end of class genome

//
//	//Calls special constructor that creates a Genome of 3 possible types:
//	//0 - Fully linked, no hidden nodes
//	//1 - Fully linked, one hidden node splitting each link
//	//2 - Fully connected with a hidden layer 
//	//num_hidden is only used in type 2
//	//Saves to file "auto_genome"
//	Genome *new_Genome_auto(int num_in,int num_out,int num_hidden,int type, const char *filename);
//
//	void print_Genome_tofile(Genome *g,const char *filename);
//
//} // namespace NEAT
//
//#endif
//
