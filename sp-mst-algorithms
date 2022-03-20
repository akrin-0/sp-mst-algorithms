/** @author Muhammet Fatih Ulu
 *  @brief A Java program for finding shortest path and minimum spanning tree of given graphs.
 */

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class sp-mst-algorithms {
	
	/**
	 * 1st part.
	 * Dijkstra's algorithm to find the shortest path of
	 * given directed graph, starting vertex and ending 
	 * vertex.
	 */
	   
	// Initialization of constants and path 
	private static final int INT_MAX = Integer.MAX_VALUE;
	private static final int NO_ROOT = -1;
	private static final int UNKNOWN = -1;
	private static String pathStr = "";
	
	// Edge class for 1st part.
	static class EdgeSP {
		int src;
		int dest;
		int distance;
		
		// Constructor.
		public EdgeSP(int src, int dest, int distance) {
			this.src = src;
			this.dest = dest;
			this.distance = distance;
		}
	}
	
	// Main class where all methods are included for 1st part.
	static class GraphSP {
		// Declare # of vertices and adjacency array list array.
		int vertices;
		ArrayList<EdgeSP>[] adjList;
		
		// Declare and initialize array lists. 
		public GraphSP(int vertices) {
			this.vertices = vertices;	
			adjList = new ArrayList[vertices];
			for (int i = 0; i < vertices; i++) {
				adjList[i] = new ArrayList<>();
			}
		}	

		// A method to add edges into array lists.
		public void addEdgeSP(int src, int dest, int weight) {
			EdgeSP edgeSp = new EdgeSP(src, dest, weight);
			adjList[src].add(edgeSp);
		}
	
		// The method to implement Dijkstra's algorithm for shortest path.
		private String[] shortestPathFinder(int startVertex, int endVertex, int vertices) {
			
			// Declarations and initializations.
			
			// Index will return the parent of index'th vertex.
			int[] parentsArray = new int[vertices];
			// Index will return shortest distance from starting vertex to index'th vertex.
			int[] shortestDist = new int[vertices];
			// Index will return true if index is visited.
			boolean[] visited = new boolean[vertices];
			
			for (int i = 0; i < vertices; i++) {
				parentsArray[i] = UNKNOWN;
				shortestDist[i] = INT_MAX;
				visited[i] = false;
			}
			
			// Shortest path from starting vertex to itself is 0.
			shortestDist[startVertex] = 0;
			// Starting vertex is root.
			parentsArray[0] = NO_ROOT;
			
			// Loop through all vertices except the starting one to find the shortest path.
			for (int i = 1; i < vertices; i++) {
	       
				// Initialize next vertex and shortest distance to that.
				int currentVertex = UNKNOWN;
				int shortestDistance = INT_MAX;
				
				// Loop through all vertices to find the current vertex to be processed.
				for (int index = 0; index < vertices; index++)
					// The vertex should be unvisited and there should be a path between vertex to another.
					if (!visited[index] && (shortestDist[index] < shortestDistance)) {
						currentVertex = index;
						shortestDistance = shortestDist[index];
					}
	                      
				// If there is no vertex to be found for index 'i', then make index 'i' visited and go next.
				if (currentVertex == UNKNOWN) {
					visited[i] = true;
					continue;
				} 
				else {			
					// Modify the current vertex as visited in the path.
					visited[currentVertex] = true;
										
					// Loop through all vertices that is connected to current vertex.
					for (int index = 0; index < adjList[currentVertex].size(); index++) {
						// Initialize edge distance from current vertex to processed index'th vertex.
						int prevEdgeDist = adjList[currentVertex].get(index).distance;
						// Value below will be used to check if it has shortest distance.
						int possibleShortestDist = shortestDistance + prevEdgeDist;
		                
						// Check if there is a shorter distance, if so update.
		                if (shortestDist[adjList[currentVertex].get(index).dest] > possibleShortestDist) {
		                	// Set current vertex as parent of looping index.
		                   	parentsArray[adjList[currentVertex].get(index).dest] = currentVertex; 
		                   	// Update to shorter path.
		            	   	shortestDist[adjList[currentVertex].get(index).dest] = possibleShortestDist;
		                }
					}
				}
			}
			
			/**
			 *  Store jumping vertex and start from 
			 *  the end. Loop and change it to its 
			 *  parent till Mecnun's city is found.
			 *  
			 *  Path string is built in reverse order,
			 *  after the loop there is a reverse again.
			 */
			int jumpingVertex = endVertex;
			while (true) {
				// If there is unknown parent of a vertex, then there is no path; so return.
				if (jumpingVertex == UNKNOWN) {
					String[] resultedStr = {"-1", "-1"};
					return resultedStr;					
				}
					
				// City "c1" is assigned to 0th index, so indexes are incremented by one. 
				String vertexStr = String.valueOf(jumpingVertex + 1);
				
				// Use StringBuilder class to reverse the number, then add 'c' and a space to that.
				StringBuilder vertexStrBld = new StringBuilder(vertexStr);
				vertexStrBld.reverse();
				pathStr += vertexStrBld + "c ";
				
				// Change the jumping vertex to check if Mecnun's city is found.
				jumpingVertex = parentsArray[jumpingVertex];
				if (jumpingVertex == 0) {
					vertexStr = String.valueOf(jumpingVertex + 1);
					vertexStrBld = new StringBuilder(vertexStr);
					vertexStrBld.reverse();
					pathStr += vertexStrBld + "c";				
					break; 
				}
			}
			// Create new StringBuilder object to reverse entire string.
			StringBuilder pathStrBld = new StringBuilder(pathStr);
			pathStrBld.reverse();
			// Create String array to carry two distinct results into main method.
			String[] resultedStr = {pathStrBld.toString(), String.valueOf(shortestDist[endVertex])};
			return resultedStr;
		}
	}
   
	
	/** 
	 * 2nd part.
	 * Kruskal's algorithm to find MST of a given 
	 * graph. It's efficient to go through edges.
	 */
   	 
	// A class for Edges.
	class EdgeMST implements Comparable<EdgeMST> {
		int src;
        int dest;
        int weight;
        /**
         *  To sort edges in Edge class array to their 
         *  weight. Less weight --> more prioritized.
         */
        public int compareTo(EdgeMST other) {
            return this.weight - other.weight;
        }
    };
 
    // Declaration of # of vertices, edges and array of edge objects.
    int verticesNum;
    int edgesNum; 
    EdgeMST[] edgesMst;
    
    // Method to create graph, initialize # of vertices, edges and Edge collection.
    void GraphMST(int vertices, int edges) {
        verticesNum = vertices;
        edgesNum = edges;
        edgesMst = new EdgeMST[edgesNum];
        for (int i = 0; i < edges; i++)
            edgesMst[i] = new EdgeMST();
    }
 
    // A class to store root and height info.
    class Root {
        int root;
        int height;
    };
    
    // A method which returns the root of given node.
    int selfRoot(Root rootsArray[], int node) {
    	// If given node is the root itself, do nothing.
    	if (rootsArray[node].root != node) 
    		// If not, find the root recursively.
    		rootsArray[node].root = selfRoot(rootsArray, rootsArray[node].root);
    	return rootsArray[node].root;
    }
 
    // A method to connect two trees (can be vertices) by using Root array.
    void connectTwo(Root rootsArray[], int node1, int node2) {
    	// Find roots of given nodes.
        int root1 = selfRoot(rootsArray, node1);
        int root2 = selfRoot(rootsArray, node2);
 
        // Check their heights. Working mechanism of heights are like in binomial queues.
        if (rootsArray[root1].height > rootsArray[root2].height)
        	rootsArray[root2].root = root1;
        else if (rootsArray[root1].height < rootsArray[root2].height)
        	rootsArray[root1].root = root2;
        else { 
        	rootsArray[root2].root = root1;
        	rootsArray[root1].height++; 
        }
    }
          
    // The method to implement Kruskal's algorithm for minimum spanning tree.
    String minSpanningTree() {
       
        // Sort the edges implemented from input to start with the shortest edge.
        Arrays.sort(edgesMst);

        // Declare and initialize the MST path.
    	EdgeMST path[] = new EdgeMST[verticesNum];
    	for (int i = 0; i < verticesNum; i++)
            path[i] = new EdgeMST();
 
    	/**
    	 *  Declare and initialize rootsArray of each 
    	 *  vertex. Initially root's of vertices are 
    	 *  themselves, and heights of them are 0.
    	 */
        Root rootsArray[] = new Root[verticesNum];
        for (int i = 0; i < verticesNum; i++) {
        	rootsArray[i] = new Root();
         	rootsArray[i].root = i;
        	rootsArray[i].height = 0;
        }
 
        // Initialize index of edge, # of edges visited and total payment made.
        int edgeIndex = 0; 
        int edgeVisited = 0;
        int minimumCost = 0;
 
        // Loop till visiting all edges.
        while (edgeVisited < (verticesNum - 1)) { // # of edges is one less than # of vertices in MST.
        	
        	/**
        	 *  If incremented index of edge is greater 
        	 *  than or equal to length of total possible 
        	 *  edges in MST, then there is no path.
        	 */        	
        	if (edgeIndex >= edgesMst.length)
        		return "-2";
        	
        	/**
        	 *  Start with the minimum cost edge and connect
        	 *  related vertices. Increase as loop continues.
        	 */
        	EdgeMST currentEdge = edgesMst[edgeIndex++];
        	
        	// Store processed edge's source's and destination's root.
            int root1 = selfRoot(rootsArray, currentEdge.src);
            int root2 = selfRoot(rootsArray, currentEdge.dest);
 
            /**
             * If edge's source's root isn't edge's destination's 
             * root, then connect them, add processed edge to
             * path and increase # of visited edges.
             */
            if (root1 != root2) { // Meaning that they are not connected.
              	connectTwo(rootsArray, root1, root2);
                path[edgeVisited++] = currentEdge;
            }
        }
        
        // Find minimum total spent by looping through the visited edges.
        for (int i = 0; i < edgeVisited; i++) {
        	minimumCost += path[i].weight;
        }
        return String.valueOf(2 * minimumCost); // Multiplication is for both Leyla's and Mecnun's travel.
    }
      
    
   
   // Main method.
   public static void main (String[] args) throws IOException {
	   	long start1 = System.nanoTime();
	    // Read all lines in input.
		ArrayList<String> allLines = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(args[0]));
			String line = reader.readLine();
			while (line != null) {
				allLines.add(line);
				line = reader.readLine();
			}				
			reader.close();
		} catch (FileNotFoundException err) {
			System.out.println("An error occured while reading the file.");
			err.printStackTrace();
		}
		
		// If the last line only consist of spaces, then size is one less than, else equal to return of size() method.
		int sizeOfLines;
		if (allLines.get(allLines.size() - 1).replaceAll("\\s+$", "").equals("")) // 
			sizeOfLines = allLines.size() - 1; // the last line is empty
		else
			sizeOfLines = allLines.size();
			
		int lineNumForPaths = sizeOfLines - 3; // The first three line of input file is fixed.
		int lineNumForPathsC = 0; // For the # of lines that is between Mecnun's and Leyla's cities.
		for (int i = 0; i < sizeOfLines; i++) { // To find the # of cities starting with c and d.
				String[] wordsArr = allLines.get(i).split(" ");
				String word = wordsArr[0];
				if (word.charAt(0) == 'c') // If first word's first char is c, then increase.
					lineNumForPathsC++;
		}
		
		/**
		 *  In third line (Mecnun's and Leyla's cities are given),
		 *  we have a line starting with 'c' which we don't want to count.
		 */  
		lineNumForPathsC--;
		
		int lineNumForPathsD = lineNumForPaths - lineNumForPathsC; // For the total lines that is starting with d.

		// Initialize time limit and total # of vertices.
		int timeLimit= Integer.parseInt(allLines.get(0));
		int verticesNum = Integer.parseInt(allLines.get(1));
		
		
		String[] wordsArr = allLines.get(2).split(" "); // In the third line, we have Mecnun's and Leyla's city.
		String mecnunCityStr = wordsArr[0]; // First word is Mecnun's city.
		String leylaCityStr = wordsArr[1];
		int mecnunCity = Integer.parseInt(mecnunCityStr.substring(1)) - 1; // 1 decrement is for e.g., to make "c1" 0th city.
		int leylaCity = Integer.parseInt(leylaCityStr.substring(1)) - 1;

		/**
		 * There is a possibility that a city can be reached 
		 * by other cities, however any other city cannot be 
		 * reached by that city; therefore it might be wrong
		 * to assign cities regarding their related # of lines.
		 * Then  we should count all cities in input by using
		 * set not to miss a city.
		 */
		
		// Creation of set of c cities.
		
		// Cities before Leyla's city.
		Set<String> set = new HashSet<String>();
		for (int i = 0; i < leylaCity; i ++) {
			String[] wordsArrSet = allLines.get(i + 3).split(" ");
			for (int j = 0; j < wordsArrSet.length; j++) {
				if (j == 0)
					set.add(wordsArrSet[j]);
				j++;
				if (j < wordsArrSet.length)
					set.add(wordsArrSet[j]);
			}			
		}
		
		// Cities after Leyla's city.
		for (int i = 0; i < (lineNumForPathsC - leylaCity - 1); i ++) {
			String[] wordsArrSet = allLines.get(i + leylaCity + 4).split(" ");
			for (int j = 0; j < wordsArrSet.length; j++) {
				if (j == 0)
					set.add(wordsArrSet[j]);
				j++;
				if (j < wordsArrSet.length)
					set.add(wordsArrSet[j]);
			}			
		}
		
		// Leyla's city line.
		int edgesInLeylaLine = 0; // Find # of edges that Leyla's city connects to d cities.
		String[] wordsArrSetL = allLines.get(leylaCity + 3).split(" ");
		for (int j = 0; j < wordsArrSetL.length; j++) {
			if (j == 0)
				set.add(wordsArrSetL[j]);
			j++;
			if (j < wordsArrSetL.length)
				if (wordsArrSetL[j].charAt(0) == 'd')
					edgesInLeylaLine++; // This is stored to find total edge # for part two.
				else
					set.add(wordsArrSetL[j]);
		}
		
		int verticesNumC = set.size(); // Total cities starting with 'c' char, will be used in first part.
		int verticesNumD = verticesNum - verticesNumC; 
		verticesNumD++; // verticesNumD is for cities starting with 'd' char, Leyla's city is included for the second part.
		
		/*
		 *  Create graph for shortest path.
		 *  Loop till Leyla's city is found,
		 *  continue looping after that line.
		 */
		GraphSP graph1st = new GraphSP(verticesNumC);
		for (int i = 0; i < leylaCity; i++) {
			int src, dest, weight;
			String[] wordsArr2 = allLines.get(i + 3).split(" ");
			src = Integer.parseInt(wordsArr2[0].substring(1)) - 1;
			for (int j = 1; j < wordsArr2.length; j++) {
				dest = Integer.parseInt(wordsArr2[j].substring(1)) - 1;
				j++;
				weight = Integer.parseInt(wordsArr2[j]);
				graph1st.addEdgeSP(src, dest, weight);
			}
		}

		// Skip Leyla's city and continue. 
		for (int i = 0; i < (lineNumForPathsC - leylaCity - 1); i++) {
			int src, dest, weight;
			String[] wordsArr2 = allLines.get(i + leylaCity + 4).split(" ");
			src = Integer.parseInt(wordsArr2[0].substring(1)) - 1;
			for (int j = 1; j < wordsArr2.length; j++) {
				dest = Integer.parseInt(wordsArr2[j].substring(1)) - 1;
				j++;
				weight = Integer.parseInt(wordsArr2[j]);
				graph1st.addEdgeSP(src, dest, weight);
			}			
		}
		
		// Initialize necessary info about first part.
		String path = "";
		int distance = -1;
		
		// Call method to find shortest part for the first part of the project.
		String[] Q1Str = graph1st.shortestPathFinder(mecnunCity, leylaCity, verticesNumC);
		
		// Get the path if it's found.
		path = Q1Str[0];
		// Get the minimum distance to check if Leyla and Mecnun marry.
		distance = Integer.parseInt(Q1Str[1]); 		
		
		// Find # of edges among d cities (Edge # between Leyla's city and d is found before).
		int edgeNum = 0;
		for (int i = 0; i < lineNumForPathsD; i++) {
			String[] wordsArr2 = allLines.get(i + lineNumForPathsC + 3).split(" ");
			for (int j = 1; j < wordsArr2.length; j++) {
				j++;
				edgeNum++;
			}
		}
		edgeNum += edgesInLeylaLine; // Total edge # to call the graph method for minimum spanning tree.
		
		// Initialize new object of main class, then use the Graph method to build graph.
		project3main graph = new project3main();
		graph.GraphMST(verticesNumD, edgeNum);
		
		/**
		 *  Since edges are considered, store and increase the index
		 *  of edges. Update the edges with the given values.
		 */
		
		// Start building second graph with Leyla's city.
		int edgeIndex = 0;
		String[] wordsArrLeylaLine = allLines.get(leylaCity + 3).split(" ");
		int srcL = 0;
		int destL, weightL;
		for (int j = 1; j < wordsArrLeylaLine.length; j++) {
			if (wordsArrSetL[j].charAt(0) == 'd') {
				destL = Integer.parseInt(wordsArrLeylaLine[j].substring(1));
				j++;
				weightL = Integer.parseInt(wordsArrLeylaLine[j]);
				graph.edgesMst[edgeIndex].src = srcL;
				graph.edgesMst[edgeIndex].dest = destL;
				graph.edgesMst[edgeIndex].weight = weightL;
				edgeIndex++;
			}
		}

		// Continue with cities starting with char 'd';
		for (int i = 0; i < lineNumForPathsD; i++) {
			int src;
			int dest;
			int weight;
			String[] wordsArr2 = allLines.get(i + lineNumForPathsC + 3).split(" ");
			src = Integer.parseInt(wordsArr2[0].substring(1));
			for (int j = 1; j < wordsArr2.length; j++) {
				dest = Integer.parseInt(wordsArr2[j].substring(1));
				j++;
				if (j >= wordsArr2.length)
					continue;
				weight = Integer.parseInt(wordsArr2[j]);
				graph.edgesMst[edgeIndex].src = src;
				graph.edgesMst[edgeIndex].dest = dest;
				graph.edgesMst[edgeIndex].weight = weight;
				edgeIndex++;
			}
		}	
		
		
		/*
		 *  Call method to find minimum spanning tree for the 
		 *  second part of the project. If all D cities are not
		 *  traveled, then minSpanningTree() method returns "-2".
		 */
		String totalSpent = graph.minSpanningTree();
		
		// Initialize answers.
		String Q1 = "";
		String Q2 = "";
		
		// If there is no path, then both answers are "-1".
		if (path.equals("-1")) {
			Q1 = "-1";
			Q2 = "-1";
		} 
		// If there is path but time limit exceeded, then they don't marry. 
		else if (distance > timeLimit) {
			Q1 = path;
			Q2 = "-1";
		} 
		// If there is path and duration is less than time limit, then they spend a honeymoon.
		else {
			Q1 = path;
			Q2 = totalSpent;
		}
		

		// Create output file and check if that file name exists.
		String outputPath = args[1];
		File output = new File(outputPath);
		try {					
			if (output.createNewFile())
				System.out.println("File is created.");
			else {
				System.out.println("File already exists.");
				return;
			}
			FileWriter writer = new FileWriter(outputPath);
			writer.write(Q1 + "\n");
			writer.write(Q2);	
			writer.close();
		} catch (Exception err) {
			System.out.println("An error occured while writing to the file.");
			err.getStackTrace();
		}
		System.out.println("Milliseconds passed: " + (System.nanoTime() - start1)/1000000);
   }
}
