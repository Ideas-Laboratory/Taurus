

#ifndef MULTILEVEL_H
#define MULTILEVEL_H


#include "solver.h"
#include "vector"
#include "set"
#include <ctime>
#include <random>
#include "SubGraph.h"
//#include "Multipole.h"
using namespace std;

#define epsilon 0.1
#define POS_SMALL_DOUBLE 1e-300
#define POS_BIG_DOUBLE   1e+300
//#define pair <float,float> DPoint
#include "RepulsiveForce.h"

#if defined(__WIN32__) || defined(_WIN32) || defined(__NT__)
#define OGDF_SYSTEM_WINDOWS
#endif
constexpr double pi = 3.14159265358979323846;
class Multilevel {
public:
    void delete_node(int del_node,int &last_selectable_index_of_S_node,vector<int> &S_node,vector<int> &position_in_node_set);
    void create_multilevel_representations(subgraph g, vector<subgraph> &sg, int &max_level, int seed = 42,
                                           int min_Graph_size = 50, int randomTries = 20);
    int get_random_node_common(int rand_index, int &last_trie_index,vector<int> &S_node,vector<int> &position_in_node_set);
    int get_random_node_with_lowest_star_mass(int rand_tries,int &last_selectable_index_of_S_node,vector<int> &mass_of_star,vector<int> &S_node,vector<int> &position_in_node_set);
    bool edgenumbersum_of_all_levels_is_linear(vector<subgraph> &sg, int act_level, int &bad_edgenr_counter);

    void init_multilevel_values(vector<subgraph> &sg,int act_level);

    void partition_galaxy_into_solar_systems(vector<subgraph> &sg, int act_level);

    void collaps_solar_systems(vector<subgraph> &sg, int act_level);

    void create_suns_and_planets(vector<subgraph> &sg, int act_level);

    void create_moon_nodes_and_pm_nodes(vector<subgraph> &sg, int act_level);

    void calculate_mass_of_collapsed_nodes(vector<subgraph> &sg,vector<double> &new_edgelength, int act_level);

    void
    create_edges_edgedistances_and_lambda_Lists(vector<subgraph> &sg, int act_level);

    void
    delete_parallel_edges_and_update_edgelength(vector<subgraph> &sg, vector<double> &new_edgelength, int act_level);

    void create_initial_placement(vector<subgraph> &sg,int act_level);

    void find_initial_placement_for_level(int level, vector<subgraph> &sg);

    void set_initial_positions_of_sun_nodes(int level, vector<subgraph> &sg);

    void set_initial_positions_of_planet_and_moon_nodes(int level, vector<subgraph> &sg, vector<int> &pm_nodes);

    void set_initial_positions_of_pm_nodes(int level, vector<subgraph> &sg, vector<int> &pm_nodes);

    void create_all_placement_sectors(int level, vector<subgraph> &sg);



    DPoint calculate_position(DPoint P, DPoint Q, float dist_P, float dist_Q);

    DPoint get_waggled_inbetween_position(DPoint s, DPoint t, double lambda);

    float norm(DPoint P, DPoint Q);

    DPoint create_random_pos(DPoint center, float radius, double angle_1,
                             double angle_2);

    DPoint get_barycenter_position(vector<DPoint> &L);

//void update_boxlength_and_cornercoordinate(vector<subgraph> &sg,int level);
//void call_FORCE_CALCULATION_step(vector<subgraph> &sg,int act_level,int max_level);
    int get_max_mult_iter(int act_level, int max_level, int node_nr);

    void set_average_ideal_edgelength(vector<subgraph> &sg, int act_level);

    void setSeed(int val);

    bool f_near_machine_precision(double distance, DPoint &force);

    double random_precision_number(double shift);

    int randomNumber(int low, int high);

    double randomDouble(double low, double high);

    float angle(DPoint a, DPoint q, DPoint r);
    float angle(float m_x,float m_y,float q_m_x,float q_m_y, float r_m_x,float r_m_y);
};
class FMMMLayout{
public:
    FMMMLayout();

    // destructor
    virtual ~FMMMLayout() { }
    virtual void call(graph &GA) ;
    void call(
            graph &GA,   //graph and layout
            const vector<double> &edgeLength); //factor for desired edge length
    double getCpuTime() {
        return time_total;
    }
    bool useHighLevelOptions() const { return m_useHighLevelOptions; }
    void useHighLevelOptions(bool uho) { m_useHighLevelOptions = uho; }
    void setSingleLevel(bool b) {m_singleLevel = b;}
    double unitEdgeLength() const { return m_unitEdgeLength; }
    void unitEdgeLength(double x) {m_unitEdgeLength = (( x > 0.0) ? x : 1);}
    bool newInitialPlacement() const { return m_newInitialPlacement; }
    void newInitialPlacement(bool nip) { m_newInitialPlacement = nip; }

    void randSeed(int p) { m_randSeed = ((0<=p) ? p : 1);}
    int randSeed() const {return m_randSeed;}
    int maxIntPosExponent() const { return m_maxIntPosExponent; }

    //! Sets the option maxIntPosExponent to \p e.
    void maxIntPosExponent(int e) {
        m_maxIntPosExponent = (((e >= 31)&&(e<=51))? e : 31);
    }
    double pageRatio() const { return m_pageRatio; }

    //! Sets the option pageRatio to \p r.
    void pageRatio(double r) {m_pageRatio = (( r > 0) ? r : 1);}
    int stepsForRotatingComponents() const { return m_stepsForRotatingComponents; }

    //! Sets the option stepsForRotatingComponents to \p n.
    void stepsForRotatingComponents(int n) {
        m_stepsForRotatingComponents = ((0<=n) ? n : 0);
    }
    double minDistCC() const { return m_minDistCC; }

    //! Sets the  minimal distance between connected components to \p x.
    void minDistCC(double x) { m_minDistCC = (( x > 0) ? x : 1);}
    int minGraphSize() const { return m_minGraphSize; }

    //! Sets the option minGraphSize to \p n.
    void minGraphSize(int n) { m_minGraphSize = ((n >= 2)? n : 2);}

    int randomTries() const { return m_randomTries; }

    //! Sets the option randomTries to \p n.
    void randomTries(int n) {m_randomTries = ((n>=1)? n: 1);}
    int maxIterFactor() const { return m_maxIterFactor; }
    double springStrength() const { return m_springStrength; }

    //! Sets the strength of the springs to \p x.
    void springStrength(double x) { m_springStrength  = ((x > 0)? x : 1);}

    //! Returns the strength of the repulsive forces.
    double repForcesStrength() const { return m_repForcesStrength; }

    //! Sets the strength of the repulsive forces to \p x.
    void repForcesStrength(double x) { m_repForcesStrength =((x > 0)? x : 1);}

    //! Sets the option maxIterFactor to \p f.
    void maxIterFactor(int f) { m_maxIterFactor = ((f>=1) ? f : 1 ); }
    double threshold() const { return m_threshold; }

    //! Sets the threshold for the stop criterion to \p x.
    void threshold(double x) {m_threshold = ((x > 0) ? x : 0.1);}
    int fixedIterations() const { return m_fixedIterations; }

    //! Sets the fixed number of iterations for the stop criterion to \p n.
    void fixedIterations(int n) { m_fixedIterations = ((n >= 1) ? n : 1);}

    //! Returns the scaling factor for the forces.
    double forceScalingFactor() const { return m_forceScalingFactor; }

    //! Sets the scaling factor for the forces to \p f.
    void forceScalingFactor(double f) { m_forceScalingFactor = ((f > 0) ? f : 1);}

    //! Returns the current setting of option coolTemperature.
    /**
     * If set to true, forces are scaled by coolValue()^(actual iteration) *
     * forceScalingFactor(); otherwise forces are scaled by forceScalingFactor().
     */
    bool coolTemperature() const { return m_coolTemperature; }

    //! Sets the option coolTemperature to \p b.
    void coolTemperature(bool b) { m_coolTemperature = b; }

    //! Returns the current setting of option coolValue.
    /**
     * This option defines the value by which forces are decreased
     * if coolTemperature is true.
     */
    double coolValue() const { return m_coolValue; }

    //! Sets the option coolValue to \p x.
    void coolValue(double x) { m_coolValue = (((x >0 )&&(x<=1) )? x : 0.99);}
    bool resizeDrawing() const { return m_resizeDrawing; }

    //! Sets the option resizeDrawing to \p b.
    void resizeDrawing(bool b) { m_resizeDrawing = b; }

    //! Returns the current setting of option resizingScalar.
    /**
     * This option defines a parameter to scale the drawing if
     * resizeDrawing() is true.
     */
    double resizingScalar() const { return m_resizingScalar; }

    //! Sets the option resizingScalar to \p s.
    void resizingScalar(double s) { m_resizingScalar = ((s > 0) ? s : 1);}

    //! Returns the number of iterations for fine tuning.
    int fineTuningIterations() const { return m_fineTuningIterations; }

    //! Sets the number of iterations for fine tuning to \p n.
    void fineTuningIterations(int n) { m_fineTuningIterations =((n >= 0) ? n : 0);}

    //! Returns the curent setting of option fineTuneScalar.
    /**
     * This option defines a parameter for scaling the forces in the
     * fine-tuning iterations.
     */
    double fineTuneScalar() const { return m_fineTuneScalar; }

    //! Sets the option fineTuneScalar to \p s
    void fineTuneScalar(double s) { m_fineTuneScalar = ((s >= 0) ? s : 1);}

    //! Returns the current setting of option adjustPostRepStrengthDynamically.
    /**
     * If set to true, the strength of the repulsive force field is calculated
     * dynamically by a formula depending on the number of nodes; otherwise the
     * strength are scaled by PostSpringStrength and PostStrengthOfRepForces.
     */
    bool adjustPostRepStrengthDynamically() const {
        return m_adjustPostRepStrengthDynamically;
    }

    //! Sets the option adjustPostRepStrengthDynamically to \p b.
    void adjustPostRepStrengthDynamically(bool b) {
        m_adjustPostRepStrengthDynamically = b;
    }

    //! Returns the strength of the springs in the postprocessing step.
    double postSpringStrength() const { return m_postSpringStrength; }

    //! Sets the strength of the springs in the postprocessing step to \p x.
    void postSpringStrength(double x) { m_postSpringStrength  = ((x > 0)? x : 1);}

    //! Returns the strength of the repulsive forces in the postprocessing step.
    double postStrengthOfRepForces() const { return m_postStrengthOfRepForces; }

    //! Sets the strength of the repulsive forces in the postprocessing step to \p x.
    void postStrengthOfRepForces(double x) {
        m_postStrengthOfRepForces = ((x > 0)? x : 1);
    }
    int  frGridQuotient() const {return m_frGridQuotient;}

    //! Sets the option frGridQuotient to \p p.
    void frGridQuotient(int p) { m_frGridQuotient = ((0<=p) ? p : 2);}
    //! Returns the current setting of option nmParticlesInLeaves.
    /**
     * Defines the maximal number of particles that are contained in
     * a leaf of the reduced bucket quadtree.
     */
    int nmParticlesInLeaves() const { return m_NMParticlesInLeaves; }

    //! Sets the option nmParticlesInLeaves to \p n.
    void nmParticlesInLeaves(int n) { m_NMParticlesInLeaves = ((n>= 1)? n : 1);}

    //! Returns the precision \a p for the <i>p</i>-term multipole expansions.
    int nmPrecision() const { return m_NMPrecision; }

    //! Sets the precision for the multipole expansions to \p p.
    void nmPrecision(int p) { m_NMPrecision  = ((p >= 1 ) ? p : 1);}



private:
    bool                  m_useHighLevelOptions; //!< The option for using high-level options.
    //FMMMOptions::PageFormatType m_pageFormat; //!< The option for the page format.
    double                m_unitEdgeLength; //!< The unit edge length.
    bool                  m_newInitialPlacement; //!< The option for new initial placement.
   // FMMMOptions::QualityVsSpeed m_qualityVersusSpeed; //!< The option for quality-vs-speed trade-off.

    //low level options
    //general options
    int                   m_randSeed; //!< The random seed.
   // FMMMOptions::EdgeLengthMeasurement m_edgeLengthMeasurement; //!< The option for edge length measurement.
   // FMMMOptions::AllowedPositions m_allowedPositions; //!< The option for allowed positions.
    int                   m_maxIntPosExponent; //!< The option for the used	exponent.

    //options for divide et impera step
    double                m_pageRatio; //!< The desired page ratio.
    int                   m_stepsForRotatingComponents; //!< The number of rotations.
   // FMMMOptions::TipOver  m_tipOverCCs; //!< Option for tip-over of connected components.
    double                m_minDistCC; //!< The separation between connected components.
   // FMMMOptions::PreSort  m_presortCCs; //!< The option for presorting connected components.

    //options for multilevel step
    bool                  m_singleLevel; //!< Option for pure single level.
    int                   m_minGraphSize; //!< The option for minimal graph size.
   // FMMMOptions::GalaxyChoice m_galaxyChoice; //!< The selection of galaxy nodes.
    int                   m_randomTries; //!< The number of random tries.

    //! The option for how to change MaxIterations.
    //! If maxIterChange != micConstant, the iterations are decreased
    //! depending on the level, starting from
    //! ((maxIterFactor()-1) * fixedIterations())
  //  FMMMOptions::MaxIterChange m_maxIterChange;

    int                   m_maxIterFactor; //!< The factor used for decreasing MaxIterations.
    //FMMMOptions::InitialPlacementMult m_initialPlacementMult; //!< The option for creating initial placement.

    //options for force calculation step
   // FMMMOptions::ForceModel m_forceModel; //!< The used force model.
    double                m_springStrength; //!< The strengths of springs.
    double                m_repForcesStrength; //!< The strength of repulsive forces.
    //FMMMOptions::RepulsiveForcesMethod m_repulsiveForcesCalculation; //!< Option for how to calculate repulsive forces.
    //FMMMOptions::StopCriterion m_stopCriterion; //!< The stop criterion.
    double                m_threshold; //!< The threshold for the stop criterion.
    int                   m_fixedIterations; //!< The fixed number of iterations for the stop criterion.
    double                m_forceScalingFactor; //!< The scaling factor for the forces.
    bool                  m_coolTemperature; //!< The option for how to scale forces.
    double                m_coolValue; //!< The value by which forces are decreased.
    //FMMMOptions::InitialPlacementForces m_initialPlacementForces; //!< The option for how the initial placement is done.

    //options for postprocessing step
    bool                  m_resizeDrawing; //!< The option for resizing the drawing.
    double                m_resizingScalar; //!< Parameter for resizing the drawing.
    int                   m_fineTuningIterations; //!< The number of iterations for fine tuning.
    double                m_fineTuneScalar; //!< Parameter for scaling forces during fine tuning.
    bool                  m_adjustPostRepStrengthDynamically; //!< The option adjustPostRepStrengthDynamically.
    double                m_postSpringStrength; //!< The strength of springs during postprocessing.
    double                m_postStrengthOfRepForces; //!< The strength of repulsive forces during postprocessing.

    //options for repulsive force approximation methods
    int                   m_frGridQuotient; //!< The grid quotient.
    //FMMMOptions::ReducedTreeConstruction m_NMTreeConstruction; //!< The option for how to construct reduced bucket quadtree.
    //FMMMOptions::SmallestCellFinding m_NMSmallCell; //!< The option for how to calculate smallest quadtratic cells.
    int                   m_NMParticlesInLeaves; //!< The maximal number of particles in a leaf.
    int                   m_NMPrecision; //!< The precision for multipole expansions.

    //other variables
    double max_integer_position; //!< The maximum value for an integer position.
    double cool_factor; //!< Needed for scaling the forces if coolTemperature is true.
    double average_ideal_edgelength; //!< Measured from center to center.
    double boxlength; //!< Holds the length of the quadratic comput. box.
    int number_of_components; //!< The number of components of the graph.
    DPoint down_left_corner; //!< Holds down left corner of the comput. box.
    vector<double> radius; //!< Holds the radius of the surrounding circle for each node.
    double time_total; //!< The runtime (=CPU-time) of the algorithm in seconds.
    FruchtermanReingold FR;
    Multilevel mul;
    void call_DIVIDE_ET_IMPERA_step(
            graph& G
);

    //! Calls the multilevel step for subGraph \p G.
    void call_MULTILEVEL_step_for_subGraph(
            subgraph& SG);

    //! Returns true iff stopCriterion() is not met
    bool running(int iter, int max_mult_iter, double actforcevectorlength);


    void call_FORCE_CALCULATION_step (
            vector<subgraph> &sg,
            int act_level,
            int max_level);

    //! Calls the postprocessing step.
    void call_POSTPROCESSING_step(
            vector<subgraph> &sg,
            int act_level,
            vector<DPoint>& F,
            vector<DPoint>& F_attr,
            vector<DPoint>& F_rep,
            vector<DPoint>& last_node_movement);


    void initialize_all_options();

    void update_low_level_options_due_to_high_level_options_settings();

    //! Imports for each node \a v of \p G its width, height and position (given from \p GA) in \p A.
    void import_NodeAttributes(
            const graph& G
);

    //! Imports for each edge e of G its desired length given via edgeLength.
    void import_EdgeAttributes (
            const graph& G,
            const vector<double>& edgeLength
);

    //! Sets the individual ideal edge length for each edge \a e.
    void init_ind_ideal_edgelength(
            const graph& G,
            int act_level);

    //! The radii of the surrounding circles of the bounding boxes are computed.
    void set_radii(const graph& G,int act_level);

    //! Exports for each node \a v in \p G_reduced the position of the original_node in \p GA.
    void export_NodeAttributes(
            graph& G_reduced,

            graph& GA);

    //! Creates a simple and loopfree copy of \p G and stores the corresponding node / edge attributes.
    /**
     * The corresponding node / edge attributes are stored in \p A_reduced and
     * \p E_reduced; the links to the copy_node and original node are stored in \p A,
     * \p A_reduced, too.
     */
    void make_simple_loopfree(
            const graph& G,
            int act_level,
            graph& G_reduced);

    //! Deletes parallel edges of \p G_reduced.
    /**
     * Saves for each set of parallel edges one representative edge in \p S and
     * saves in \p new_edgelength the new edge length of this edge in \p G_reduced.
     */
 /*   void delete_parallel_edges(
            const graph& G,
            EdgeArray<EdgeAttributes>& E,
            graph& G_reduced,
            List<edge>& S,
            vector<double>& new_edgelength);*/

    //! Sets for each edge \a e of \a G_reduced in \p S its edgelength to \p new_edgelength[\a e].
    /**
     * Also stores this information in \p E_reduced.
     */
    /*void update_edgelength(
            vector<edge>& S,
            vector<double>& new_edgelength,
            EdgeArray<EdgeAttributes>& E_reduced);*/

    //! Returns the value for the strength of the repulsive forces.
    /**
     * Used in the postprocessing step; depending on \p n = G.numberOfNodes().
     */
    double get_post_rep_force_strength(int n) {
        return min(0.2,400.0/double(n));
    }

    void adjust_positions( vector<subgraph>& sg,
                          int act_level);


    //! @}
    //! \name Functions for divide et impera step
    //! @{

    //! Constructs the list of connected components of G.
    /**
     * Also constructs the corresponding lists with the node / edge attributes
     * (containing a pointer to the original node in \p G for each node in a subgraph).
     */
    void create_maximum_connected_subGraphs(
            graph& G,
            graph G_sub[],
            vector<int>& component);

    //! The drawings of the subgraphs are packed.
    /**
     * This is done such that the subgraphs do not overlap and fit into a small
     * box with the desired aspect ratio.
     */
    void pack_subGraph_drawings(
            graph G_sub);

    //! The bounding rectangles of all connected componenents of \a G are calculated and stored in \p R.
   /* void  calculate_bounding_rectangles_of_components(
            List<Rectangle>& R,
            Graph  G_sub[],
            NodeArray<NodeAttributes> A_sub[]);

    //! The bounding rectangle of the componenet_index-th. component of G is returned.
    Rectangle calculate_bounding_rectangle(
            Graph& G,
            NodeArray<NodeAttributes>& A,
            int componenet_index);*/

    /**
     * If number_of_components > 1, the subgraphs \p G_sub are rotated and skipped to
     * find bounding rectangles with minimum area. The information is saved in \p R and
     * the node positions in \p A_sub are updated. If number_of_components == 1 a rotation
     * with minimal aspect ratio is found instead.
     */
/*
    void rotate_components_and_calculate_bounding_rectangles(
            List<Rectangle>&R,
            Graph G_sub[],
            NodeArray<NodeAttributes> A_sub[]);
*/

    /**
     * Returns the area (aspect ratio area) of a rectangle with width w and height h
     * if comp_nr > 1 ( comp_nr == 1).
     */
    double calculate_area(double width, double height, int comp_nr) {
        double scaling = 1.0;
        if (comp_nr == 1) {  //calculate aspect ratio area of the rectangle
            if( height == 0.0 ){ fprintf(stderr,"height=0 error\n");}
            double ratio = width / height;
            if (ratio < pageRatio()) { //scale width
                if( ratio == 0.0 ){fprintf(stderr,"ratio=0 error\n");}
                scaling = pageRatio() / ratio;
            } else { //scale height
                if( pageRatio() == 0.0 ){fprintf(stderr,"pageRatio=0 error\n");}
                scaling = ratio / pageRatio();
            }
        }
        return width * height * scaling;
    }


    void calculate_forces(
            vector<subgraph> &sg,
            int act_level,
            vector<DPoint>& F,
            vector<DPoint>& F_attr,
            vector<DPoint>& F_rep,
            vector<DPoint>& last_node_movement,
            int iter,
            int fine_tuning_step);

    //! The length of the computational box in the first iteration is set (down left corner is at (0,0).
    void init_boxlength_and_cornercoordinate(vector<subgraph> &sg,int level);

    //! The initial placements of the nodes are created by using initialPlacementForces().
    void create_initial_placement (vector<subgraph>& sg,
                                   int act_level);



    //! Sets all entries of \p F to (0,0).
    void  init_F (vector<subgraph>& sg,
                  int act_level, vector<DPoint>& F);


    //! Make initializations for the data structures that are used in the choosen class for rep. force calculation.
    void make_initialisations_for_rep_calc_classes(
            vector<subgraph>& sg,
            int act_level/*,
		NodeArray<NodeAttributes> &A,
		NodeArray<DPoint>& F_rep*/);

    //! Calculates repulsive forces for each node.
    void calculate_repulsive_forces(
            vector<subgraph>& sg,
            int act_level,
            vector<DPoint>& F_rep);
   /* {
               // NM.calculate_repulsive_forces(G,A,F_rep);

    }*/


    //! Deallocates dynamically allocated memory of the choosen rep. calculation class.
    void deallocate_memory_for_rep_calc_classes()
    {

           // NM.deallocate_memory();
    }

    //! Calculates attractive forces for each node.
    void calculate_attractive_forces(
            vector<subgraph>& sg,
            int act_level,
            vector<DPoint>& F_attr);

    //! Returns the attractive force scalar.
    double f_attr_scalar (double d,double ind_ideal_edge_length);

    //! Add attractive and repulsive forces for each node.
    void add_attr_rep_forces(
            vector<subgraph>& sg,
            int act_level,
            vector<DPoint>& F_attr,
            vector<DPoint>& F_rep,
            vector<DPoint>& F,
            int iter,
            int fine_tuning_step);

    //! Move the nodes.
    void move_nodes(vector<subgraph>& sg,
                    int act_level,vector<DPoint>& F);

    //! Computes a new tight computational square-box.
    /**
     * (Guaranteeing, that all midpoints are inside the square.)
     */
    void update_boxlength_and_cornercoordinate(vector<subgraph> &sg,int level);

    //! Describes the max. radius of a move in one time step, depending on the number of iterations.
    double max_radius(int iter) {
        return (iter == 1) ? boxlength/1000 : boxlength/5;
    }

    //! The average_ideal_edgelength for all edges is computed.
    void set_average_ideal_edgelength(vector<subgraph> &sg,int act_level);

    /**
     * Calculates the average force on each node in the actual iteration, which is
     * needed if StopCriterion is scThreshold() or scFixedIterationsOrThreshold().
     */
    double get_average_forcevector_length (vector<subgraph>& sg,
                                           int act_level, vector<DPoint>& F);

    /**
     * Depending on the direction of \p last_node_movement[\a v], the length of the next
     * displacement of node \a v is restricted.
     */
    void prevent_oscillations(
            vector<subgraph>& sg,
            int act_level,
            vector<DPoint>& F,
            vector<DPoint>&
            last_node_movement,
            int iter);

    //! \p last_node_movement is initialized to \p F (used after first iteration).
    void init_last_node_movement(
            vector<subgraph>& sg,
            int act_level,
            vector<DPoint>& F,
            vector<DPoint>& last_node_movement);

    /**
     * If resizeDrawing is true, the drawing is adapted to the ideal average
     * edge length by shrinking respectively expanding the drawing area.
     */
    void adapt_drawing_to_ideal_average_edgelength(
            vector<subgraph>& sg,
            int act_level);

    /**
     * The force is restricted to have values within the comp. box (needed for
     * exception handling, if the force is too large for further calculations).
     */
    void restrict_force_to_comp_box(DPoint& force) {
        double x_min = down_left_corner.first;
        double x_max = down_left_corner.first+boxlength;
        double y_min = down_left_corner.second;
        double y_max = down_left_corner.second+boxlength;
        if (force.first < x_min )
            force.first = x_min;
        else if (force.first > x_max )
            force.first = x_max;
        if (force.second < y_min )
            force.second = y_min;
        else if (force.second > y_max )
            force.second = y_max;
    }

};

#endif
