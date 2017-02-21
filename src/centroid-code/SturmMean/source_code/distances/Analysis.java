package distances;

import polyAlg.PolyMain;
import static polyAlg.PolyMain.getGeodesic;
import static polyAlg.PolyMain.calcGeoDist;
import polyAlg.Tools;
import distanceAlg1.PhyloTree;
import distanceAlg1.Geodesic;
import distanceAlg1.PhyloTreeEdge;
import distanceAlg1.EdgeAttribute;
import distanceAlg1.Bipartition;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import centroid.*;


public class Analysis {
	
	public static int verbose = 0;

	
	
	
	
	/**  Returns the average geodesic distance for all the geodesics in geos.
	 * 
	 * @param geos
	 * @return
	 */
	public static double averageGeoDist(Geodesic[][] geos) {
		int numTrees = geos.length;
		double average = 0;
		int numDists = 0;
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				average = average + (geos[i][j]).getDist();
				numDists++;
			}
		}
		return average / numDists;
	}
	
	/** Returns the average distance of the trees from the origin.
	 *  i.e. the average lengths/norms of the trees
	 * @param trees
	 * @return
	 */
	public static double avgDistanceFromOrigin(PhyloTree [] trees) {
		int numTrees = trees.length;
		double avg = 0;
		
		for(int i = 0; i < trees.length; i++) {
			avg = avg + trees[i].getDistanceFromOrigin();
		}
		return avg/numTrees;
	}
	
	
	/** Given a set of trees, removes all edges common to these trees, and returns the new trees.
	 *  TODO:  Finish method
	 *  
	 * @param inTrees
	 * @return
	 */
	public static PhyloTree[] removeCommonSplits(PhyloTree[] inTrees) {
		PhyloTree[] outTrees = new PhyloTree[inTrees.length];
		
		/* Fill in  */  // TODO
		Vector<Bipartition> commonSplits = getCommonSplits(inTrees);
		
		if (commonSplits == null) {
			return inTrees;
		}
		
		if (commonSplits.size() == 0) {
			return inTrees;
		}
		
		for (int i = 0; i < inTrees.length;i++) {
			outTrees[i] = new PhyloTree(inTrees[i]);
			outTrees[i].removeSplits(commonSplits);
		}
		
		return outTrees;
	}	
		
	
	/** Returns a vector of all edges common to the trees in trees.  Returns null if no common edges.
	 * 
	 * @param trees
	 * @return
	 */
		
	// XXX: This method assumes that there are no 0 length edges.  However, this is currently (June 2011)
	// not compatible with the rest of the code (but shouldn't cause a problem in most cases.)
	// TODO:  assumes PhyloTree contains the same number of trees as its length
	public static Vector<Bipartition> getCommonSplits(PhyloTree[] trees) {
		Vector<Bipartition> commonSplits = new Vector<Bipartition>();
				
		if (trees == null) {
			System.err.println("No trees for finding common edges: returning null.");
			return null;
		}
		
		if (trees[0] == null) {
			System.err.println("No trees for finding common edges:  returning null");
			return null;
		}
		
		commonSplits = trees[0].getSplits();
		
		if (trees.length ==1) {
			System.out.println("One tree:  Returning all its edges as common");
			return commonSplits; 
		}
		
		for (int i = 1; i < trees.length; i++) {
			commonSplits.retainAll(trees[i].getSplits());
		}
		return commonSplits;
	}
		
		
		
	/**  Do numRep bootstraps of finding the centroid of a set of trees sampled randomly with replacement from trees.
	 * 	 Do at most numIter iterations of Sturm's algorithm to determine the centroid at each iteration.
	 *   If the algorithm hasn't converged by then, centroid will be null;
	 * @param trees
	 * @param numRep
	 * @param numIter
	 * @return
	 */
	public static PhyloTree[] bootstrap(PhyloTree[] trees, int numRep, long numIter) {
		PhyloTree[] centroids = new PhyloTree[numRep];
		int numTrees = trees.length;
		
		PhyloTree[] sampledTrees = new PhyloTree[numTrees];
		
		int counter = 0;  //  remove all statements with // next to them.
		
		for (int i = 0; i< numRep; i++) {
			sampledTrees = sampleWithReplacement(trees);
			
			double epsilon = CentroidMain.getEpsilon(sampledTrees,5000);
			
			centroids[i] = CentroidMain.getCentroidViaRandomCauchy(sampledTrees, numIter,5, epsilon, null,0);	
			
			double var = CentroidMain.variance(centroids[i],sampledTrees);  //
			
			PhyloTree starTree = getStarTree(trees); //  remove all statements with // next to them.
			
			double starVar = CentroidMain.variance(starTree, sampledTrees);  //
			
			if (var <= starVar) {  //
				System.out.println("------------------------------------  STAR NOT THE MEAN ---------------- ");  //
				System.out.println("Star is " + starTree.getNewick(true));  //
				System.out.println("Centroid is " + centroids[i].getNewick(true));   //
				System.out.println("Star tree variance is " + starVar);    //
				System.out.println("Centroid variance is " + var);    //
				
				double dist = calcGeoDist(centroids[i], starTree); //
				System.out.println("Distance to star tree is " + dist);  //
				counter++;
			}    //
		}
		
		System.out.println(counter + " bootstraps where star didn't have a lower variance");
		
		return centroids;
	}
	
	
	
	/**  Returns the proportion of geodesics that are the cone path.
	 * 
	 * @param geos
	 * @return
	 */
	public static double proportionConePath(Geodesic[][] geos) {
		int numTrees = geos.length;
		int numConePath = 0;
		int numDists = 0;
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				if (geos[i][j].getRS().size() == 1) {
					numConePath++;
				}
				numDists++;
			}
		}
		return (double) numConePath / numDists;
	}
	

	public static void geodesicStats(PhyloTree[] trees) {
		int numTrees = trees.length;
		double averageDist = 0;
		int numDists = 0;
		int numLeaves = trees[0].getLeaf2NumMap().size();
		
		double[][] dist = new double[numTrees][numTrees];
		
		/*  Initialize an array to hold the number of geodesics with the given number of ratios in its ratio sequence. */
		int[] numRSOfEachSize = new int[numLeaves];
		Arrays.fill(numRSOfEachSize,0);
		
		/* Initialize an array to hold the number of geodesics with the given number of common edges. */
		int[] numCommonEdgesOverGeos = new int[numLeaves];
		Arrays.fill(numCommonEdgesOverGeos,0);
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				Geodesic geo = PolyMain.getGeodesic(trees[i], trees[j], null);
				dist[i][j] = geo.getDist();
				
				numRSOfEachSize[geo.getRS().size()]++;
				numCommonEdgesOverGeos[geo.numCommonEdges()]++;
				
				averageDist = averageDist + geo.getDist();
				numDists++;
			}
		}
		
		double[] starDist = distancesToStarTree(trees);
		
				
		printDistancesToFile(dist,starDist, "montpellier_analysis0to19_2.txt");
		
		System.out.println("Average length of geodesic is " + averageDist/numDists);
		
		for (int i = 1; i < numTrees; i++) {
			System.out.println("# of geodesics with " + i + " support pairs: " + numRSOfEachSize[i]);
		}
		
		System.out.println();
		
		for (int i = 0; i < numTrees; i++) {
			System.out.println("# of geodesics with " + i + " common edges: " + numCommonEdgesOverGeos[i]);
		}
	}
	
	
	
	/**  Returns the star tree with averaged leaf edge lengths.
	 * 
	 * @param trees
	 * @return
	 */
	public static PhyloTree getStarTree(PhyloTree[] trees) {
		if (trees == null) {
			System.err.println("Warning:  tree set null when getting star tree; returning null");
			return null;
		}
		
		if (trees[0] == null) {
			System.err.println("Warning:  getting star tree for empty set of trees; returning null");
			return null;
		}
		
		int numLeaves = trees[0].getLeaf2NumMap().size();
		int numTrees = trees.length;
		
		Vector<PhyloTreeEdge> starEdges = new Vector<PhyloTreeEdge>();
		
		PhyloTree star = new PhyloTree(starEdges,trees[0].getLeaf2NumMap(), trees[0].isRooted()); 
		
		EdgeAttribute[] starLeafAttribs = new EdgeAttribute[numLeaves];
		
		
		for (int leaf= 0; leaf < numLeaves; leaf++) {
			starLeafAttribs[leaf] = trees[0].getLeafEdgeAttribs()[leaf].clone();
			for (int tree = 1; tree < numTrees; tree++) {
				if (starLeafAttribs[leaf] == null) {
					System.out.println("starLeafAttribs is null for leaf " + leaf);
				}
				else if (trees[tree] == null) {
					System.out.println("tree " + tree + " is null");
				}
				else if (trees[tree].getLeafEdgeAttribs()[leaf] == null) {
					System.out.println("trees[tree].getLeafEdgeAttribs()[leaf] is null for tree " + tree + " and leaf " + leaf);
				}
				starLeafAttribs[leaf] = EdgeAttribute.add(starLeafAttribs[leaf], trees[tree].getLeafEdgeAttribs()[leaf].clone());
			}
			starLeafAttribs[leaf].scaleBy(1.0/numTrees);
		}
		star.setLeafEdgeAttribs(starLeafAttribs);
		
		return star;
	}
	
	
	public static double[] distancesToStarTree(PhyloTree[] trees) {
		int numTrees = trees.length;
		
		PhyloTree star = getStarTree(trees);
		
		double[] dist = new double[numTrees];
		
		for (int i = 0; i < numTrees; i++) {
			dist[i] = calcGeoDist(star,trees[i]);
		}
		return dist;
	}
	
	
	
	public static void printDistancesToFile(double[][] dist, double[] specialDist, String fileName) {
		int numTrees = dist.length;
		
		PrintWriter outputStream = null;
  
        /* Outputs the distances in a matrix, with 0 entries on the diagonals, and the rows labelled from 1 to numTrees. */
        try {
        	outputStream = new PrintWriter(new FileWriter(fileName));
 
    		for (int i = 0; i < numTrees; i++) {
    			outputStream.print(i);
    			
    			if (i != 0) {
    				outputStream.print("\t" + specialDist[i-1]);
    			}
    			
    			for (int j = 1; j< i; j++) {
    				outputStream.print("\t" + dist[j-1][i-1]);
    			}
    			outputStream.println("\t0");
    		}
    		if (outputStream != null) {
    			outputStream.close();
    		}
        } catch (FileNotFoundException e) {
        	System.out.println("Error opening or writing to " +fileName + ": "+ e.getMessage());
        	System.exit(1);
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + fileName + ": " + e.getMessage());
        	System.exit(1);
        }
	}
	
	/** Write the lengths of the trees in trees to the file outfile.
	 * 
	 * @param trees
	 * @param outfile
	 */
	public static void writeLengths(PhyloTree[] trees, String outfile) {		
		
		PrintWriter outputStream = null;
  
        /* Outputs the distances in a matrix, with 0 entries on the diagonals, and the rows labelled from 1 to numTrees. */
        try {
        	outputStream = new PrintWriter(new FileWriter(outfile));
 
    		for (int i = 0; i < trees.length; i++) {
    			
    			outputStream.println(trees[i].getDistanceFromOrigin());
    			
    		}
    		if (outputStream != null) {
    			outputStream.close();
    		}
        } catch (FileNotFoundException e) {
        	System.out.println("Error opening or writing to " +outfile + ": "+ e.getMessage());
        	System.exit(1);
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + outfile + ": " + e.getMessage());
        	System.exit(1);
        }
	}
	
	
	/** Returns an array of the trees of the same size as trees, with the entries drawn uniformly from trees with replacement.
	 * 
	 * 
	 * @param trees
	 * @return
	 */
	public static PhyloTree[] sampleWithReplacement(PhyloTree[] trees) {
		int numTrees = trees.length;
		PhyloTree[] sampledTrees = new PhyloTree[numTrees];
		
		XORShiftRandom rng = new XORShiftRandom();
		
		
		for (int i = 0; i < numTrees; i++) {
			sampledTrees[i] = trees[rng.nextInt(numTrees)].clone();
		}
		
		return sampledTrees;		
	}
	
	
	public static void main(String[] args) {
		String outfile = "output.txt";
		String treeFile = "";
		String otherTreeFile = null;
		Boolean rooted = true;
		String algorithm = "";
		PhyloTree[] otherTrees = null;
		double epsilon = -1;
		
		/* Parse command line arguments */

		if (args.length < 1) {
			System.out.println("Error: Missing input file name"); displayHelp(); System.exit(1);
		}
		treeFile = args[args.length-1];
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); System.exit(1); }
			// specify algorithm
			else if (args[i].equals("-a")) {
				if (i < args.length -2) { algorithm = args[i+1]; i++; }
				else { System.err.println("Error: algorithm not specified"); displayHelp(); System.exit(1); }
			}
			// epsilon, if desired
			else if (args[i].equals("-e")) {
				if (i < args.length -2) { epsilon = Double.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: value for epsilon not specified"); displayHelp(); System.exit(1); }
			}
			// other tree file, if desired
			else if (args[i].equals("-f")) {
				if (i < args.length -2) { otherTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of other tree file not specified"); displayHelp(); System.exit(1); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); System.exit(1); }
			}	
				
			// all other arguments.  Note we can have -vn
			else { 
				for (int j = 1; j<args[i].length(); j++) {
					switch(args[i].charAt(j)) {						
					// display help
					case 'h': displayHelp(); System.exit(0); break;
					case 'u': rooted = false; break;
					
					default: System.out.println("Illegal command line option.\n"); displayHelp(); System.exit(1); break;
					} // end switch
				} // end j loop (arguments without parameter)
			} // end parsing an individual argument
		}  // end for i (looping through arguments)
		
		/* Read in the trees  */
		
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		
		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}
		
		if (algorithm.equals("compute_SSD")) {
			// otherTreeFile should contain a single tree
			// treeFile contains the list of trees to compute sum of square distance to
			if (otherTreeFile == null) {
				System.out.println("Error:  second tree file not specified");
				System.exit(1);
			}
			
			
			double ssd = CentroidMain.sumOfSquareDist(otherTrees[0], trees);
			System.out.println("Sum of squares distance is " + ssd);
			
		}
		else if (algorithm.equals("gtp_twofiles")) {
			// otherTreeFile contains one list of trees to compare against the trees in treeFile
			if (otherTreeFile == null) {
				System.out.println("Error:  second tree file not specified");
				System.exit(1);
			}
			
			computeAllGeodesicsBtwLists(trees,otherTrees,outfile);
		}
		
		else if (algorithm.equals("verify_treefile")) {
			PolyMain.readInTreesFromFile(treeFile, rooted);
			System.out.println("Done reading in trees from file " + treeFile);
			System.exit(0);
		}
		
		//  Compute projection of first tree in otherTreeFile onto the geodesic between
		// the first two trees in treeFile
		else if (algorithm.equals("project")) {
			projectTreeToGeo(trees,otherTrees,epsilon,outfile);
			System.exit(0);
		}
		
		// Writes the lengths of the trees in treefile to outfile, 
		// one per line
		else if (algorithm.equals("lengths")) {
			writeLengths(trees, outfile);
			System.exit(0);
		}
		
		//  Nicely prints out the splits that differ in the trees in trees from the trees in otherTrees
		else if (algorithm.equals("diff_splits")) {
			printIncompatibleSplits(trees,otherTrees,outfile, false);
			System.exit(0);
		}
		
	//  Nicely prints out the splits that differ between tree i in trees and tree i in otherTrees
		else if (algorithm.equals("diff_splits_paired")) {
			printIncompatibleSplits(trees,otherTrees,outfile, true);
			System.exit(0);
		}
		
		else {
			System.out.println("Error:  no algorithm specified.\n");
			System.exit(1);
		}
		
		System.exit(0);
		
	}
	
	/** Help message (ie. which arguments can be used, etc.)
	 * 
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("java -jar analysis.jar [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -a <algorithm> \t specifies what to compute");
		System.out.println("\t -e <epsilon> \t specifies epsilon, if needed.");
		System.out.println("\t -f <otherTreeFile> \t reads in an additional tree file");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>.  Default is output.txt");
		System.out.println("\t -u \t trees are unrooted. Default is trees are rooted.");
		System.out.println("\n");
		System.out.println("Algorithms are:");
		System.out.println("\t compute_SSD \t computes the sum of square distance from the tree in otherTreeFile to the trees in treefile");
		System.out.println("\t gtp_twofiles \t computes the geodesic distance between all trees in treefile and all trees in otherTreeFile");
		System.out.println("\t verify_treefile \t verifies that all trees in the treefile can be read in without errors");
		System.out.println("\t project \t projects the first tree in otherTreeFile onto the geodesic between the first two trees in treefile");
		System.out.println("\t lengths \t computes the lengths (norm of edge vectors) of all trees in treefile");
		System.out.println("\t diff_splits \t nicely prints out (to the output file) the splits that differ between trees in treefile and trees in otherTreeFile");
		System.out.println("\t diff_splits_paired \t nicely prints out (to the output file) the splits the differ between tree i in treefile and tree i in otherTreeFile");
	}
	
	/**  Computes all geodesics between the two lists.
	 *   Prints the results to file.
	 * 
	 * @param trees
	 * @param otherTrees
	 * @param rooted
	 * @param outfile
	 */
	public static void computeAllGeodesicsBtwLists(PhyloTree[] trees,PhyloTree[] otherTrees,String outFileName) {
		int numTrees = trees.length;
		int numOtherTrees = otherTrees.length;
		
//		Geodesic[][] geos = new Geodesic[numTrees][numOtherTrees]; 
//		double [][] dists = new double[numTrees][numOtherTrees];
		           
//		for (int i = 0; i < numTrees; i++) {
//			for (int j = 0; j < numOtherTrees; j++) {
//				geos[i][j] = PolyMain.getGeodesic(trees[i], otherTrees[j], "geo_" + i + "_" + j);
//				dists[i][j] = geos[i][j].getDist();
//			}
//		}
	    
	    // print distances to file
	    PrintWriter outputStream = null;
	  
	    // Outputs the distances in a column, with the first two columns being the trees numbers and the third
	    // number the geodesic distance between those trees
	    try {
	       	outputStream = new PrintWriter(new FileWriter(outFileName));
	 
	    	for (int i = 0; i < numTrees ; i++) {
	    		for (int j = 0; j< numOtherTrees; j++) {
	    			double dist = PolyMain.getGeodesic(trees[i], otherTrees[j], "geo_" + i + "_" + j).getDist();
//					outputStream.println(i + "\t" + j + "\t" + Tools.round(geos[i][j].getDist(), 6));
	    			outputStream.println(i + "\t" + j + "\t" + Tools.round(dist, 6));
				}
				outputStream.println();
			}
			if (outputStream != null) {
	            outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFileName + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	    	System.out.println("Error opening or writing to " + outFileName + ": " + e.getMessage());
	    	System.exit(1);
	    }
	}
	
	/** Wrapper for code that projects one tree onto the geodesic bewteen two others.
	 * Computes the projection of otherTrees[0] onto the geodesic between trees[0] and trees[1].
	 * Writes the projected tree to outFileName.
	 * 
	 * @param trees
	 * @param otherTrees
	 * @param outFileName
	 */
	public static void projectTreeToGeo(PhyloTree[] trees, PhyloTree[] otherTrees, double epsilon, String outFileName) {
		if ((trees.length < 2) || (otherTrees.length < 1)) {
			System.out.println("Error: not enough trees in input files");
			System.exit(1);
		}
		
		if (epsilon <= 0) {
			epsilon = 0.05;
			System.out.println("Using epsilon = 0.05 when projecting tree onto geodesic");
		}
		
		PhyloTree treeToProject = otherTrees[0];
		PhyloTree t1 = trees[0];
		PhyloTree t2 = trees[1];
		
		PhyloTree projectedTree = treeToProject.projectToGeo(getGeodesic(t1,t2,null),epsilon); 
		
		 // print projected tree to file
	    PrintWriter outputStream = null;
	  
	    try {
	       	outputStream = new PrintWriter(new FileWriter(outFileName));
	 
			outputStream.println(projectedTree.getNewick(true));
			
			if (outputStream != null) {
	            outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFileName + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	    	System.out.println("Error opening or writing to " + outFileName + ": " + e.getMessage());
	    	System.exit(1);
	    }
	}

	
	/** Prints nicely to a file info about which splits in the trees in otherTree in are incompatible
	 *  with the trees in goodTrees (considered individually).
	 * If goodTrees contains only one tree (i.e. an original tree), then this prints 
	 * information about which splits in the trees in otherTrees are incompatible with it.
	 * 
	 */
	public static void printIncompatibleSplits(PhyloTree[] goodTrees, PhyloTree[] otherTrees, String outFilename, Boolean paired) {
		String output = "";
		Vector<String> leaf2NumMap = goodTrees[0].getLeaf2NumMap();
		int startIndex, endIndex;
		
		for (int i = 0; i < goodTrees.length; i++) {
			output = output + "------------------------------------\n";
			output = output + "------------------------------------\n";
			output = output + "Tree " + i + " in input file:\n";
			
			
			if (paired) {
				startIndex = i;
				endIndex = i + 1;
			}
			else {
				startIndex = 0;
				endIndex = otherTrees.length;
			}
			
			for (int j = startIndex; j < endIndex; j++) {
				output = output + "\n";
				
				Vector<PhyloTreeEdge> incompEdges = goodTrees[i].getEdgesIncompatibleWith(otherTrees[j]);
				
//				System.out.println("otherTree edges incomp with good tree: ");
				
				if (incompEdges.size() >0) {
					output = output + "\tOther tree " + j + " is missing splits:\n";
					// display edges in verbose output, in a numbered list
					for (int k = 0; k < incompEdges.size(); k++) {
						output = output + "\t" + k + ". " + incompEdges.get(k).toStringReroot(leaf2NumMap, "Opossum") + "\n";
//						System.out.println(incompEdges.get(k).toStringVerbose(leaf2NumMap));
					}
					
				}
				
				incompEdges = otherTrees[j].getEdgesIncompatibleWith(goodTrees[i]);
				
				if (incompEdges.size() > 0) {
					output = output + "\n\thas extra splits:\n";
					
					// display edges in verbose output, in a numbered list
					for (int k = 0; k < incompEdges.size(); k++) {
						output = output + "\t" + k + ". " + incompEdges.get(k).toStringVerbose(leaf2NumMap) + "\n";
					}
				}
				
			}
			output = output + "\n";
		}
		
		
		// write the output to the file
		// write verbose output to geofile, if in verbose mode
	    PrintWriter outputStream = null;
	        
	    try {
	        outputStream = new PrintWriter(new FileWriter(outFilename));
	        	
	        outputStream.println(output);
	        
	    	if (outputStream != null) {
	               outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFilename + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	        System.out.println("Error opening or writing to " + outFilename + ": " + e.getMessage());
	        System.exit(1);
	    }
	}
}
