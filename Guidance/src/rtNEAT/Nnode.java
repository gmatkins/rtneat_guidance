package rtNEAT;

import java.util.Vector;

//#ifndef _NNODE_H_
//#define _NNODE_H_
//
//#include <algorithm>
//#include <vector>
//#include "neat.h"
//#include "trait.h"
//#include "link.h"
//
//namespace NEAT {
class Nnode{
//
//	enum nodetype {
//		NEURON = 0,
//		SENSOR = 1
//	};
	public enum nodetype {
		NEURON, SENSOR
	}
	
//
//	enum nodeplace {
//		HIDDEN = 0,
//		INPUT = 1,
//		OUTPUT = 2,
//		BIAS = 3
//	};
	public enum nodeplace{
		HIDDEN, INPUT, OUTPUT, BIAS
	}
//
//	enum functype {
//		SIGMOID = 0
//	};
	public enum functype {
		SIGMOID
	}
//
//	class Link;
//	
//	class Network;
//
//	// ----------------------------------------------------------------------- 
//	// A NODE is either a NEURON or a SENSOR.  
//	//   - If it's a sensor, it can be loaded with a value for output
//	//   - If it's a neuron, it has a list of its incoming input signals (List<Link> is used) 
//	// Use an activation count to avoid flushing
//	class NNode {
//
//		friend class Network;
//		friend class Genome;
//
//	protected:
//
//		int activation_count;  // keeps track of which activation the node is currently in
	protected int activation_count;
//		double last_activation; // Holds the previous step's activation for recurrency
	protected double last_activation;
//		double last_activation2; // Holds the activation BEFORE the prevous step's
	protected double last_activation2;
//		// This is necessary for a special recurrent case when the innode
//		// of a recurrent link is one time step ahead of the outnode.
//		// The innode then needs to send from TWO time steps ago
//
//		Trait *nodetrait; // Points to a trait of parameters
	protected Trait nodetrait;
//
//		int trait_id;  // identify the trait derived by this node
	protected int trait_id;
//
//		NNode *dup;       // Used for Genome duplication
	protected Nnode dup;
//
//		NNode *analogue;  // Used for Gene decoding
	protected Nnode analogue;
//
//		bool override; // The NNode cannot compute its own output- something is overriding it
	protected boolean override;
//
//		double override_value; // Contains the activation value that will override this node's activation
	protected double override_value;
//
//		// Pointer to the Sensor corresponding to this Body.
//		//Sensor* mySensor;
//
//	public:
//		bool frozen; // When frozen, cannot be mutated (meaning its trait pointer is fixed)
	public boolean frozen;
//
//		functype ftype; // type is either SIGMOID ..or others that can be added
	public functype ftype;
//		nodetype type; // type is either NEURON or SENSOR
	public nodetype type;
//
//		double activesum;  // The incoming activity before being processed
	public double activesum;
//		double activation; // The total activation entering the NNode 
	public double activation;
//		bool active_flag;  // To make sure outputs are active
	public boolean active_flag;
//
//		// NOT USED IN NEAT - covered by "activation" above
//		double output;  // Output of the NNode- the value in the NNode 
	public double output;
//
//		// ************ LEARNING PARAMETERS *********** 
//		// The following parameters are for use in    
//		//   neurons that learn through habituation,
//		//   sensitization, or Hebbian-type processes  
//
//		double params[NEAT::num_trait_params];
	public double[] params = new double[Neat.NUM_TRAIT_PARAMS];
//
//		std::vector<Link*> incoming; // A list of pointers to incoming weighted signals from other nodes
	public Vector<Link> incoming;
//		std::vector<Link*> outgoing;  // A list of pointers to links carrying this node's signal
	public Vector<Link> outgoing;
//
//		// These members are used for graphing with GTK+/GDK
//		std::vector<double> rowlevels;  // Depths from output where this node appears
	public Vector<Double> rowlevels;
//		int row;  // Final row decided upon for drawing this NNode in
	public int row;
//		int ypos;
	public int ypos;
//		int xpos;
	public int xpos;
//
//		int node_id;  // A node can be given an identification number for saving in files
	public int node_id;
//
//		nodeplace gen_node_label;  // Used for genetic marking of nodes
	public nodeplace gen_node_label;
//
//		NNode(nodetype ntype,int nodeid);
	public Nnode(nodetype ntype,int nodeid) {
		active_flag=false;
		activesum=0;
		activation=0;
		output=0;
		last_activation=0;
		last_activation2=0;
		type=ntype; //NEURON or SENSOR type
		activation_count=0; //Inactive upon creation
		node_id=nodeid;
		ftype=functype.SIGMOID;
		gen_node_label= nodeplace.HIDDEN;
		frozen=false;
		trait_id=1;
		override=false;
	}
//
//		NNode(nodetype ntype,int nodeid, nodeplace placement);
	public Nnode(nodetype ntype,int nodeid, nodeplace placement) {
		active_flag=false;
		activesum=0;
		activation=0;
		output=0;
		last_activation=0;
		last_activation2=0;
		type=ntype; //NEURON or SENSOR type
		activation_count=0; //Inactive upon creation
		node_id=nodeid;
		ftype=functype.SIGMOID;
		gen_node_label=placement;
		frozen=false;
		trait_id=1;
		override=false;
	}
//
//		// Construct a NNode off another NNode for genome purposes
//		NNode(NNode *n,Trait *t);
	public Nnode(Nnode n, Trait t) {
		active_flag=false;
		activation=0;
		output=0;
		last_activation=0;
		last_activation2=0;
		type=n.type; //NEURON or SENSOR type
		activation_count=0; //Inactive upon creation
		node_id=n.node_id;
		ftype=functype.SIGMOID;
		gen_node_label=n.gen_node_label;
		nodetrait=t;
		frozen=false;
		if (t != null)
			trait_id=t.trait_id;
		else trait_id=1;
		override=false;
	}
//
//		// Construct the node out of a file specification using given list of traits
//		NNode (const char *argline, std::vector<Trait*> &traits);
	Nnode ( String argline, Vector<Trait> traits) {
		int traitnum;
		//std::vector<Trait*>::iterator curtrait;
		
		activesum=0;

	    //std::stringstream ss(argline);
		//char curword[128];
		//char delimiters[] = " \n";
		//int curwordnum = 0;

		//Get the node parameters
		//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
		//node_id = atoi(curword);
		//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
		//traitnum = atoi(curword);
		//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
		//type = (nodetype)atoi(curword);
		//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
		//gen_node_label = (nodeplace)atoi(curword);

	    int nodety, nodepl;
	    //ss >> node_id >> traitnum >> nodety >> nodepl;
	    String[] splitArgline = argline.split(" ");
	    node_id = Integer.parseInt(splitArgline[0]);
	    traitnum = Integer.parseInt(splitArgline[1]);
	    nodety = Integer.parseInt(splitArgline[2]);
	    nodepl = Integer.parseInt(splitArgline[3]);
	    
	    //type = (nodetype)nodety;
	    type = nodetype.values()[nodety];
	    gen_node_label = nodeplace.values()[nodepl];

		// Get the Sensor Identifier and Parameter String
		// mySensor = SensorRegistry::getSensor(id, param);
		frozen=false;  //TODO: Maybe change

		//Get a pointer to the trait this node points to
		if (traitnum==0) nodetrait = null;
		else {
			//curtrait=traits.begin();
			//while(((*curtrait)->trait_id)!=traitnum)
			//	++curtrait;
			//nodetrait=(*curtrait);
			//trait_id=nodetrait->trait_id;
			for (Trait curtrait : traits){
				nodetrait = curtrait;
				trait_id = nodetrait.trait_id;
			}
		}

		override=false;
	}
//
//		// Copy Constructor
//		NNode (const NNode& nnode);
	public Nnode (Nnode nnode)
	{
		active_flag = nnode.active_flag;
		activesum = nnode.activesum;
		activation = nnode.activation;
		output = nnode.output;
		last_activation = nnode.last_activation;
		last_activation2 = nnode.last_activation2;
		type = nnode.type; //NEURON or SENSOR type
		activation_count = nnode.activation_count; //Inactive upon creation
		node_id = nnode.node_id;
		ftype = nnode.ftype;
		nodetrait = nnode.nodetrait;
		gen_node_label = nnode.gen_node_label;
		dup = nnode.dup;
		analogue = nnode.dup;
		frozen = nnode.frozen;
		trait_id = nnode.trait_id;
		override = nnode.override;
	}
//
//		~NNode();
//
//		// Just return activation for step
//		double get_active_out();
	
//
//		// Return activation from PREVIOUS time step
//		double get_active_out_td();
//
//		// Returns the type of the node, NEURON or SENSOR
//		const nodetype get_type();
	public nodetype get_type() {
		return type;
	}
//
//		// Allows alteration between NEURON and SENSOR.  Returns its argument
//		nodetype set_type(nodetype);
	nodetype set_type(nodetype newtype) {
		type=newtype;
		return newtype;
	}
//
//		// If the node is a SENSOR, returns true and loads the value
//		bool sensor_load(double);
	public boolean sensor_load(double value) {
		if (type==nodetype.SENSOR) {

			//Time delay memory
			last_activation2=last_activation;
			last_activation=activation;

			activation_count++;  //Puts sensor into next time-step
			activation=value;
			return true;
		}
		else return false;
	}
//
//		// Adds a NONRECURRENT Link to a new NNode in the incoming List
//		void add_incoming(NNode*,double);
	public void add_incoming(Nnode feednode,double weight,boolean recur) {
		Link newlink=new Link(weight,feednode,this,recur);
		incoming.add(newlink);
		(feednode.outgoing).add(newlink);
	}
//
//		// Adds a Link to a new NNode in the incoming List
//		void add_incoming(NNode*,double,bool);
	public void add_incoming(NNode *feednode,double weight) {
		Link *newlink=new Link(weight,feednode,this,false);
		incoming.push_back(newlink);
		(feednode->outgoing).push_back(newlink);
	}
//
//		// Recursively deactivate backwards through the network
//		void flushback();
//
//		// Verify flushing for debugging
//		void flushback_check(std::vector<NNode*> &seenlist);
//
//		// Print the node to a file
//        void  print_to_file(std::ostream &outFile);
//	void print_to_file(std::ofstream &outFile);
//
//		// Have NNode gain its properties from the trait
//		void derive_trait(Trait *curtrait);
//
//		// Returns the gene that created the node
//		NNode *get_analogue();
//
//		// Force an output value on the node
//		void override_output(double new_output);
//
//		// Tell whether node has been overridden
//		bool overridden();
//
//		// Set activation to the override value and turn off override
//		void activate_override();  
//
//		// Writes back changes weight values into the genome
//		// (Lamarckian trasnfer of characteristics)
//		void Lamarck();
//
//		//Find the greatest depth starting from this neuron at depth d
//		int depth(int d,Network *mynet); 
//
//	};
//
//
//} // namespace NEAT
	
} // endl Nnode class

//
//#endif
//
