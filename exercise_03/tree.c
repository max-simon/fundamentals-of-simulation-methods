#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>


/* Let's define two types of structures, one for the particles, one for the tree nodes
 */

typedef struct particle_data
{
  double pos[3];
  double mass;

  double acc_tree[3];       //acceleration of particle calculated from tree
  double acc_exact[3];      //acceleration of particle calculated by exact summation
}
  particle;


typedef struct node_data
{
  double center[3];    // geom. center
  double len;

  double cm[3];         // center of mass
  double mass;

  struct node_data *suns[8];    //subnodes

  particle *p;
}
  node;


#define  MAX_POINTS               40000      /* this sets the particle number that is used */
static double opening_threshold = 0.8;      /* tree opening angle */
static double eps               = 0.001;    /* gravitational softening length */


/* Lets create, for simplicity, some static arrays which will hold our data */

#define  MAX_NODES    (5 * MAX_POINTS)

node      tree[MAX_NODES];
particle  star[MAX_POINTS];


/* this function returns a pointer to an empty tree node
 */
node *get_empty_node(void)
{
  node *no;
  static int count_nodes = 0;

  if(count_nodes < MAX_NODES)
    {
      no = &tree[count_nodes++];        //first free entry of tree
      memset(no, 0, sizeof(node));      //set node "no" to 0
    }
  else
    {
      printf("sorry, we are out of tree nodes.\n");
      exit(1);
    }

  return no;
}


/* this function determines the index (0-7) of one of the 8 subnodes in which
 * a particle falls within a given node
 */
int get_subnode_index(node *current, particle *p)
{
  int index = 0;

  if(p->pos[0] > current->center[0])
    index += 4;
  if(p->pos[1] > current->center[1])
    index += 2;
  if(p->pos[2] > current->center[2])
    index += 1;

  return index;
}


/* this function's task is to insert a new particle into a given tree node
 */
void insert_particle(node *current, particle *pnew)
{
  node *sun;
  int i, j, k, n, p_subnode, pnew_subnode;

  if(current->p)
    {
      /* The node contains a particle already.
	 Need to create a new set of 8 subnodes, and then move this particle to one of them */

      for(i = 0, n = 0; i<2; i++)
      for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    {
	      sun = get_empty_node();
	      current->suns[n++] = sun;
	      sun->len = 0.5 * current->len;
	      sun->center[0] = current->center[0] + 0.25 * (2*i-1) * current->len;
	      sun->center[1] = current->center[1] + 0.25 * (2*j-1) * current->len;
	      sun->center[2] = current->center[2] + 0.25 * (2*k-1) * current->len;
	    }

      /* determine in which subnode the old particle sits */
      p_subnode = get_subnode_index(current, current->p);

      /* now move the particle to this subnode */
      current->suns[p_subnode]->p = current->p;
      current->p = NULL;

      /* determine in which subnode the new particle sits */
      pnew_subnode = get_subnode_index(current, pnew);

      /* now try to insert the new particle there */
      insert_particle(current->suns[pnew_subnode], pnew);
    }
  else
    {
      /* check in which subnode the new particle would fall */
      pnew_subnode = get_subnode_index(current, pnew);

      /* if the corresponding subnode exists, we try to insert the particle there,
	 otherwise we know there are no subnodes in the node, so we can put the particle into the current node */
      if(current->suns[pnew_subnode])
	insert_particle(current->suns[pnew_subnode], pnew);
      else
	current->p = pnew;
    }
}


/* This function recursively calculates the multipole moments for the current node.
 * We only use monopoles here.
 */
void calc_multipole_moments(node *current)
{
  int n, j;

  if(current->suns[0])   /* do we have subnodes? */
    {
      /* yes, so let's first calculate their multipole moments */
      for(n = 0; n < 8; n++)
      calc_multipole_moments(current->suns[n]);

      /* initialize the node multipole moments to zero */
      current->mass  = 0;
      for(int j = 0; j < 3; j++)
       current->cm[j] = 0;

      /* now calculate the moment of the current node from those of its suns */

      /* Done */
      for(int n = 0; n < 8; n++) { // Loop over all subnodes
        current->mass += current->suns[n]->mass; // calculate total mass
        for(int j = 0; j < 3; j++) {
          current->cm[j] += current->suns[n]->mass * current->suns[n]->cm[j]; // calculate M * center of mass
        }
      }
      for(int j = 0; j < 3; j++) {
        current->cm[j] = current->cm[j]/current->mass; // calculate center of mass
      }

    }
  else
    {
      if(current->p)  /* do we at least have a particle? */
	{
	  /* yes, so let's copy this particle to the multipole moments of the node */

	  current->mass = current->p->mass;
	  for(j = 0; j < 3; j++)
	    current->cm[j] = current->p->pos[j];
	}
      else
	{
	  /* nothing in here at all; let's initialize the multipole moments to zero */
	  current->mass  = 0;
	  for(j = 0; j < 3; j++)
	    current->cm[j] = 0;
	}
    }
}


double get_opening_angle(node *current, double pos[3])
{
  int j;
  double r2 = 0;

  for(j = 0; j < 3; j++)
    r2 += (current->cm[j] - pos[j]) * (current->cm[j] - pos[j]);

  return current->len / (sqrt(r2) + 1.0e-35);       //"+ 1.0e-35" avoids divison by 0
}



double walk_tree(node *current, double pos[3], double acc[3], double counter)
{
  int n;
  double theta;


  if(current->mass)   /* only do something if there is mass in this branch of the tree (i.e. if it is not empty) */
    {
      theta = get_opening_angle(current, pos);

      /* if the node is seen under a small enough angle or contains a single particle,
       * we take its multipole expansion, and we're done for this branch
       */
      if(theta < opening_threshold || current->p)
	{
	  /*
	   Done
	   */
     double y[3]; // calculate y-vector
     for(int i = 0; i < 3; i++) {
       y[i] = pos[i] - current->cm[i];
     }
     double norm_y = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);

     for(int i = 0; i < 3; i++) {
       acc[i] -= current->mass * y[i]/pow((norm_y*norm_y) + (eps * eps), 3/2); // see exercise sheet for formula
     }
     counter += 1; // count, how many nodes were used

	}
      else
	{
	  /* otherwise we have to open the node and look at all daughter nodes in turn */

	  if(current->suns[0])             /* make sure that we actually have subnodes */
	    for(n=0; n<8; n++)
	      counter = walk_tree(current->suns[n], pos, acc, counter);
	}
    }
  return counter;
}



int main(int argc, char **argv)
{
  node *root;
  int i, j, N;
  double t0, t1;

  N = MAX_POINTS;

  srand48(42);   /* set a random number seed */

  /* create a random particle set, uniformly distributed in a box */
  for(i=0; i < N; i++)
    {
      star[i].mass = 1.0 / N;       // total mass = 1

      for(j=0; j<3; j++)
	star[i].pos[j] = drand48();
    }

  /* create an empty root node for the tree */
  root = get_empty_node();

  /* set the dimension and position of the root node */
  root->len = 1.0;
  for(j=0; j<3; j++)
    root->center[j] = 0.5;

  /* insert the particles into the tree */
  for(i=0; i < N; i++)
    insert_particle(root, &star[i]);

  /* calculate the multipole moments */
  calc_multipole_moments(root);


  /* set a timer */
  t0 = (double) clock();

  double counter = 0; // count, how many nodes were opened (for all particles)

  /* now calculate the accelerations with the tree */
  for(i = 0; i < N; i++)
    {
      for(j = 0; j < 3; j++) {
	     star[i].acc_tree[j] = 0;
     }

      counter = walk_tree(root, star[i].pos, star[i].acc_tree, counter);
    }

  double mean_used_nodes = counter/N;
  t1 = (double) clock();
  double time_tree = (t1 - t0) / CLOCKS_PER_SEC;
  printf("\nforce calculation with tree took:        %8g sec\n", time_tree);
  printf("\nper particle %8g nodes were used", mean_used_nodes);


  t0 = (double) clock();

  /* now calculate the accelerations with direct summation, for comparison */
  for(int i = 0; i < N; i++)
    {
	  /*
	   Done
	   */
     for(int k = 0; k < 3; k++) { // initialize to 0
       star[i].acc_exact[k] = 0;
     }
     for(int j = 0; j < N; j++) { // loop over all other particles
       if(j == i) { // if it is the same particle
         continue;
       }
       // calculate distance
       double norm_d = sqrt(pow(star[i].pos[0] - star[j].pos[0], 2) + pow(star[i].pos[1] - star[j].pos[1], 2) + pow(star[i].pos[2] - star[j].pos[2], 2));
       for(int k = 0; k < 3; k++) {
         star[i].acc_exact[k] -= star[j].mass * (star[i].pos[k] - star[j].pos[k])/pow(norm_d*norm_d + (eps * eps), 3/2); // see sheet.pdf for formula
       }
     }
    }

  t1 = (double) clock();
  double time_exact = (t1 - t0) / CLOCKS_PER_SEC;
  printf("\ncalculation with direct summation took:  %8g sec\n", time_exact);

  /* now do the calculation of the mean relative error
   */

  /* Done */
  double err_sum = 0;
  double diff[3] = {0};
  double norm_exact = 0;
  double norm_diff = 0;
  for(int i = 0; i < N; i++) {
    norm_exact = sqrt(pow(star[i].acc_exact[0], 2) + pow(star[i].acc_exact[1], 2) + pow(star[i].acc_exact[2], 2));
    for(int j = 0; j < 3; j++) {
      diff[j] = star[i].acc_exact[j] - star[i].acc_tree[j];
    }
    norm_diff = sqrt(pow(diff[0], 2) + pow(diff[1], 2) + pow(diff[2], 2));
    err_sum += norm_diff/norm_exact;
  }
  err_sum /= N;

  printf("\nAverage relative error:  %8g\n", err_sum);


  // Append results to file
  FILE * f = fopen("results.txt", "a");

  fprintf(f, "%d  %1.3f  %8g  %8g  %8g  %8g\n", N, opening_threshold, time_tree, mean_used_nodes, time_exact, err_sum);

  fclose(f);

  exit(0);
}
