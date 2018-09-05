#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#ifdef __linux__
#include <execinfo.h>
#endif
#include <unistd.h>
#include <omp.h>

/*#include <utils.h>
#include <constant.h>
#include <integral.h>
#include <particle.h>
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

/// gcc main.c -o main `pkg-config --cflags --libs gsl` -lm -Wall

/** Compilation with GSL */
/// gcc main.c -o omp_parallel -fopenmp `pkg-config --cflags --libs gsl` -lquadmath -lm -Wall


/** Compilation with the external library */
/// gcc -Wall -o main main.c `pkg-config --cflags --libs gsl` -L /home/cayley/Desktop/Thesis/Programmation/AAA_Library -Wl,-rpath=/home/cayley/Desktop/Thesis/Programmation/AAA_Library -I /home/cayley/Desktop/Thesis/Programmation/AAA_Library -lutils -lconstant -lintegral -lparticle -lm


/**  Constant CGS */
#define sigmaT 6.6524e-25
#define c 2.99792458e10
#define c2 8.9875517873681764e20
#define qe 4.8032068e-10
#define me 9.1093897e-28
#define mec2 8.1871111680068e-7
#define mp 1.6749286e-24
#define Msun 1.989e33
#define h 6.6260755e-27
#define h2 4.390487653170025e-53
#define hbar 1.05457266e-27
#define kb 1.380657e-16
#define sigmaSB 5.67051e-5
#define ath 7.5657e-15
#define Ggrav 6.67259e-8
#define Bcrit 4.414e13


/** Constant for Cosmology */
#define Hubble 69.6              /// km s^-1 Mpc^-1
#define Mpc_to_cm 3.08568e24               /// km s^-1 Mpc^-1
#define s_to_year 3.17098e-8

/** Unit conversion */
#define Hz_to_eV 4.13566553855e-15      /// Convert Hz to eV

/** Constant for kinetic theory */
#define alphaf 0.0072973525664
#define lambdac 2.4263102367e-10
#define hplank 6.6261e-27
#define cyclosynchrotron 2.7992483604657531e6          /// e/(2 pi m_e c)   x B


/**  Constant m_e c^2 */
#define me_u 1.0
#define mmu_u 206.768289933



/** Numbers */
#define sqrt2 1.4142135623730950
#define sqrt3 1.7320508075688773
#define pi    3.1415926535897932
#define pi2   9.8696044010893586
#define Z3    1.202056903


#define GRID_E 100
#define GRID_P 150



double vabs(double x){
    if(x<0.0){
        return -x;
    }
    else{
        return x;
    }

}

/***** Integral */
const double GAUSS_NODE_7[7] = { -0.949107912342758524526, -0.7415311855993944398639, -0.4058451513773971669066, 0.0, 0.4058451513773971669066, 0.7415311855993944398639, 0.949107912342758524526};
const double GAUSS_WEIGHT_7[7] = {  0.129484966168869693271, 0.279705391489276667901, 0.38183005050511894495, 0.4179591836734693877551, 0.38183005050511894495, 0.279705391489276667901, 0.129484966168869693271};
const double KRONROD_NODE_15[15] = {   -0.991455371120812639207, -0.949107912342758524526, -0.86486442335976907279, -0.7415311855993944398639, -0.5860872354676911302941, -0.4058451513773971669066, -0.2077849550078984676007, 0.0, 0.2077849550078984676007, 0.4058451513773971669066, 0.5860872354676911302941, 0.7415311855993944398639, 0.86486442335976907279, 0.949107912342758524526, 0.991455371120812639207};
const double KRONROD_WEIGHT_15[15] = { 0.0229353220105292249637, 0.063092092629978553291,  0.10479001032225018384, 0.140653259715525918745, 0.1690047266392679028266, 0.1903505780647854099133, 0.204432940075298892414, 0.209482141084727828013, 0.204432940075298892414, 0.1903505780647854099133, 0.1690047266392679028266, 0.140653259715525918745,  0.10479001032225018384, 0.063092092629978553291, 0.0229353220105292249637};
const double GAUSS_NODE_15[15] = {  -0.9879925180204854284895657,  -0.9372733924007059043077589, -0.8482065834104272162006483, -0.72441773136017004741618605, -0.57097217260853884753722674, -0.39415134707756336989720737, -0.2011940939974345223006283, 0.0, 0.2011940939974345223006283, 0.3941513470775633698972074,  0.5709721726085388475372267, 0.7244177313601700474161861, 0.8482065834104272162006483, 0.9372733924007059043077589, 0.9879925180204854284895657};
const double GAUSS_WEIGHT_15[15] = { 0.03075324199611726835462839, 0.07036604748810812470926742, 0.1071592204671719350118695, 0.13957067792615431444780479, 0.16626920581699393355320086, 0.1861610000155622110268006, 0.1984314853271115764561183, 0.2025782419255612728806202, 0.19843148532711157645611833,  0.18616100001556221102680056, 0.1662692058169939335532009, 0.1395706779261543144478047, 0.1071592204671719350118695, 0.0703660474881081247092674, 0.0307532419961172683546284};

const double KRONROD_NODE_51[51] = { -9.992621049926098341934574865403406e-01, -9.955569697904980979087849468939016e-01, -9.880357945340772476373310145774062e-01, -9.766639214595175114983153864795941e-01, -9.616149864258425124181300336601672e-01, -9.429745712289743394140111696584705e-01, -9.207471152817015617463460845463306e-01, -8.949919978782753688510420067828050e-01, -8.658470652932755954489969695883401e-01, -8.334426287608340014210211086935696e-01, -7.978737979985000594104109049943066e-01, -7.592592630373576305772828652043610e-01, -7.177664068130843881866540797732978e-01, -6.735663684734683644851206332476222e-01, -6.268100990103174127881226816245179e-01, -5.776629302412229677236898416126541e-01, -5.263252843347191825996237781580102e-01, -4.730027314457149605221821150091920e-01, -4.178853821930377488518143945945725e-01, -3.611723058093878377358217301276407e-01, -3.030895389311078301674789099803393e-01, -2.438668837209884320451903627974516e-01, -1.837189394210488920159698887595284e-01, -1.228646926107103963873598188080368e-01, -6.154448300568507888654639236679663e-02, 0.0, 6.154448300568507888654639236679663e-02, 1.228646926107103963873598188080368e-01, 1.837189394210488920159698887595284e-01, 2.438668837209884320451903627974516e-01, 3.030895389311078301674789099803393e-01, 3.611723058093878377358217301276407e-01, 4.178853821930377488518143945945725e-01, 4.730027314457149605221821150091920e-01, 5.263252843347191825996237781580102e-01, 5.776629302412229677236898416126541e-01, 6.268100990103174127881226816245179e-01, 6.735663684734683644851206332476222e-01, 7.177664068130843881866540797732978e-01, 7.592592630373576305772828652043610e-01, 7.978737979985000594104109049943066e-01, 8.334426287608340014210211086935696e-01, 8.658470652932755954489969695883401e-01, 8.949919978782753688510420067828050e-01, 9.207471152817015617463460845463306e-01, 9.429745712289743394140111696584705e-01, 9.616149864258425124181300336601672e-01, 9.766639214595175114983153864795941e-01, 9.880357945340772476373310145774062e-01, 9.955569697904980979087849468939016e-01,  9.992621049926098341934574865403406e-01};
const double KRONROD_WEIGHT_51[51] = { 1.987383892330315926507851882843410e-03, 5.561932135356713758040236901065522e-03, 9.473973386174151607207710523655324e-03, 1.323622919557167481365640584697624e-02, 1.684781770912829823151666753633632e-02, 2.043537114588283545656829223593897e-02, 2.400994560695321622009248916488108e-02, 2.747531758785173780294845551781108e-02, 3.079230016738748889110902021522859e-02, 3.400213027432933783674879522955120e-02, 3.711627148341554356033062536761988e-02, 4.008382550403238207483928446707565e-02, 4.287284502017004947689579243949516e-02, 4.550291304992178890987058475266039e-02, 4.798253713883671390639225575691475e-02, 5.027767908071567196332525943344008e-02, 5.236288580640747586436671213787271e-02, 5.425112988854549014454337045987561e-02, 5.595081122041231730824068638274735e-02, 5.743711636156783285358269393950647e-02, 5.868968002239420796197417585678776e-02, 5.972034032417405997909929193256185e-02, 6.053945537604586294536026751756543e-02, 6.112850971705304830585903041629271e-02, 6.147118987142531666154413196526418e-02, 6.158081806783293507875982424006455e-02, 6.147118987142531666154413196526418e-02, 6.112850971705304830585903041629271e-02, 6.053945537604586294536026751756543e-02, 5.972034032417405997909929193256185e-02, 5.868968002239420796197417585678776e-02, 5.743711636156783285358269393950647e-02, 5.595081122041231730824068638274735e-02, 5.425112988854549014454337045987561e-02, 5.236288580640747586436671213787271e-02, 5.027767908071567196332525943344008e-02, 4.798253713883671390639225575691475e-02, 4.550291304992178890987058475266039e-02, 4.287284502017004947689579243949516e-02, 4.008382550403238207483928446707565e-02, 3.711627148341554356033062536761988e-02, 3.400213027432933783674879522955120e-02, 3.079230016738748889110902021522859e-02, 2.747531758785173780294845551781108e-02, 2.400994560695321622009248916488108e-02, 2.043537114588283545656829223593897e-02, 1.684781770912829823151666753633632e-02, 1.323622919557167481365640584697624e-02, 9.473973386174151607207710523655324e-03, 5.561932135356713758040236901065522e-03, 1.987383892330315926507851882843410e-03};
const double GAUSS_NODE_25[25] = { -9.955569697904980979087849468939016e-01, -9.766639214595175114983153864795941e-01, -9.429745712289743394140111696584705e-01, -8.949919978782753688510420067828050e-01, -8.334426287608340014210211086935696e-01, -7.592592630373576305772828652043610e-01, -6.735663684734683644851206332476222e-01, -5.776629302412229677236898416126541e-01, -4.730027314457149605221821150091920e-01, -3.611723058093878377358217301276407e-01, -2.438668837209884320451903627974516e-01, -1.228646926107103963873598188080368e-01, 0.0, 1.228646926107103963873598188080368e-01, 2.438668837209884320451903627974516e-01, 3.611723058093878377358217301276407e-01, 4.730027314457149605221821150091920e-01, 5.776629302412229677236898416126541e-01, 6.735663684734683644851206332476222e-01, 7.592592630373576305772828652043610e-01, 8.334426287608340014210211086935696e-01, 8.949919978782753688510420067828050e-01, 9.429745712289743394140111696584705e-01, 9.766639214595175114983153864795941e-01, 9.955569697904980979087849468939016e-01};
const double GAUSS_WEIGHT_25[25] = { 1.139379850102628794790296411323477e-02, 2.635498661503213726190181529529914e-02, 4.093915670130631265562348771164595e-02, 5.490469597583519192593689154047332e-02, 6.803833381235691720718718565670797e-02, 8.014070033500101801323495966911130e-02, 9.102826198296364981149722070289165e-02, 1.005359490670506442022068903926858e-01, 1.085196244742636531160939570501166e-01, 1.148582591457116483393255458695558e-01, 1.194557635357847722281781265129010e-01, 1.222424429903100416889595189458515e-01, 1.231760537267154512039028730790501e-01, 1.222424429903100416889595189458515e-01, 1.194557635357847722281781265129010e-01, 1.148582591457116483393255458695558e-01, 1.085196244742636531160939570501166e-01, 1.005359490670506442022068903926858e-01, 9.102826198296364981149722070289165e-02, 8.014070033500101801323495966911130e-02, 6.803833381235691720718718565670797e-02, 5.490469597583519192593689154047332e-02, 4.093915670130631265562348771164595e-02, 2.635498661503213726190181529529914e-02, 1.139379850102628794790296411323477e-02};
const double GAUSS_NODE_51[51] = { -0.998909990848903495168995877273386, -0.994261260436752574621084897949263, -0.98591599173590299658388570755831, -0.973903368019323867231755486394187, -0.958267848613908194557707038316324, -0.939067544002962383435367806390905, -0.916373862309780230823571294251476, -0.890271218029527303277795370736727, -0.860856711182292371473495743716111, -0.828239763823064832854818424016562, -0.792541712099381205234410878375876, -0.753895354485375525763960025452474, -0.712444457577036644580524855400146, -0.668343221175370086864460419403989, -0.621755704600723273755042745403316, -0.5728552163513038365223947025901909, -0.521823669366185842514087784826818, -0.468850904286041063610457258811622, -0.4141339832263038779368718097446579, -0.3578764566884095097752010885196664, -0.3002876063353319395302456496444203, -0.241581666447798703846733114869262, -0.181977026957077545323998701169214, -0.121695421018888766963820420963181, -0.0609611001505787247341947068432054, 0.0, 0.0609611001505787247341947068432054, 0.1216954210188887669638204209631811, 0.1819770269570775453239987011692144, 0.2415816664477987038467331148692624, 0.3002876063353319395302456496444203, 0.3578764566884095097752010885196664, 0.4141339832263038779368718097446579, 0.4688509042860410636104572588116225, 0.521823669366185842514087784826818, 0.572855216351303836522394702590191, 0.621755704600723273755042745403316, 0.668343221175370086864460419403989, 0.712444457577036644580524855400146, 0.753895354485375525763960025452474, 0.792541712099381205234410878375876, 0.828239763823064832854818424016562, 0.860856711182292371473495743716111, 0.890271218029527303277795370736727, 0.916373862309780230823571294251476, 0.939067544002962383435367806390905, 0.958267848613908194557707038316324, 0.973903368019323867231755486394187, 0.98591599173590299658388570755831, 0.994261260436752574621084897949263, 0.998909990848903495168995877273386};
const double GAUSS_WEIGHT_51[51] = { 0.002796807171089895575544216881809, 0.0065003377832526002921093751985375, 0.01018519129782172993923759178660673, 0.0138326340064778222966884530218564, 0.017428714723401052259503646433292, 0.020959988401703210579792618401514, 0.02441330057378143427314156541689821, 0.0277757985941624771959956656322089, 0.031034971290160008454425502956547, 0.0341786932041883362362093385947675, 0.037195268923260292842908275811871, 0.040073476285496453186809115921396, 0.042802607997880086653609514244286, 0.045372511407650068748166814988438, 0.047773626240623101999995353707354, 0.04999702015005740977954885536200579, 0.052034421936697087564136447468662, 0.053878252313045561434099301696972, 0.055521652095738693016737059093624, 0.056958507720258662100077726734277, 0.058183473982592140598437877661776, 0.0591919939229615437835390077491546, 0.059980315777503252090063987996517, 0.060545506934737795138125251467754, 0.06088546484485634388119861422269621, 0.060998924841205880159797643098356, 0.060885464844856343881198614222696, 0.060545506934737795138125251467754, 0.059980315777503252090063987996517, 0.059191993922961543783539007749155, 0.0581834739825921405984378776617759, 0.056958507720258662100077726734277, 0.055521652095738693016737059093624, 0.053878252313045561434099301696972, 0.052034421936697087564136447468662, 0.0499970201500574097795488553620058, 0.0477736262406231019999953537073536, 0.045372511407650068748166814988438, 0.042802607997880086653609514244286, 0.040073476285496453186809115921396, 0.037195268923260292842908275811871, 0.034178693204188336236209338594768, 0.0310349712901600084544255029565474, 0.027775798594162477195995665632209, 0.024413300573781434273141565416898, 0.020959988401703210579792618401514, 0.017428714723401052259503646433292, 0.013832634006477822296688453021856, 0.0101851912978217299392375917866067, 0.0065003377832526002921093751985375, 0.002796807171089895575544216881809};



/**  Structure for 1D integral refinment */
typedef struct refinment_1D{
    double a;
    double b;

    double I2;
    double I1;
    double I2g;

    int refine;
    int refineI2g;

    struct refinment_1D *previous;
    struct refinment_1D *next;

}refinment_1D;


/**  Refinment method */
int rewind_str_1D(refinment_1D **As){
    while((*As)->previous != NULL){
        (*As) = (*As)->previous;
    }
    return 0;
}

int Counting_1D(refinment_1D **As){

    int i;
    i = 1;
    rewind_str_1D(As);
    while((*As)->next != NULL){
        i +=1;
        (*As) = (*As)->next;
    }
    return i;
}

int refine_grid_1D(refinment_1D **As, double (*f)(double, double *), double *param, int loga, int order){

    rewind_str_1D(As);
    //printf("After rewind\n");
    int counter = 0;
    while(counter != 1){
        if((*As)->refine == 1){

            //printf("refine_grid_1D 1\n");
            /// Refinment is needed.
            double I1 = 0.0;
            double I2 = 0.0;
            double f_eval;
            int i;

            double a, b, barc;
            a = (*As)->a;
            b = (*As)->b;

            refinment_1D *temporaire = NULL;
            temporaire = (refinment_1D*) malloc(sizeof(refinment_1D));
            //printf("refine_grid_1D 2\n");
            //getchar();

            (*temporaire).refine = 0;
            //printf("refine_grid_1D 2a\n");
            (*temporaire).refineI2g = 1;
            //printf("refine_grid_1D 2b\n");
            (*As)->refineI2g = 1;
             //printf("refine_grid_1D 2c\n");
            (*temporaire).previous = (*As);
            (*temporaire).next = (*As)->next;
            (*As)->next = temporaire;
            //printf("refine_grid_1D 2d\n");
            //getchar();
            if((*temporaire).next !=NULL){
                (*((*temporaire).next)).previous =temporaire;
            }
            //printf("refine_grid_1D 3\n");


            /// proceed to the refinment itself now:
            /// 1) for As:
            if(loga == 0){
                (*As)->a = a;
                barc = (a+b)/2.0;
                (*As)->b = barc;
            }
            else{
                (*As)->a = a;
                barc = exp((log(b)+log(a))/2.0);
                (*As)->b = barc;
            }


            I1 = 0.0;
            I2 = 0.0;
            /**  Gauss - Kronrod 15 - 7 */
            if(order == 15){
                for(i=0;i<15;i++){
                    f_eval = f((barc-a)*KRONROD_NODE_15[i]/2.0+(a+barc)/2.0,param);
                    I2+= f_eval*KRONROD_WEIGHT_15[i];
                    if(i%2 ==1){
                        I1+= f_eval*GAUSS_WEIGHT_7[(i-1)/2];
                    }
                }
            }
            /**  Gauss - Kronrod 51 - 25 */
            if(order == 51){
                for(i=0;i<51;i++){
                    f_eval = f((barc-a)*KRONROD_NODE_51[i]/2.0+(a+barc)/2.0,param);
                    I2+= f_eval*KRONROD_WEIGHT_51[i];
                    if(i%2 ==1){
                        I1+= f_eval*GAUSS_WEIGHT_25[(i-1)/2];
                    }
                }
            }

            I2 = I2*(barc-a)/2.0;
            I1 = I1*(barc-a)/2.0;
            (*As)->I1 = I1;
            (*As)->I2 = I2;

            /// for temporaire
            (*temporaire).a = barc;            /// The log or non log refinement is already taking care in the c, computed above
            (*temporaire).b = b;

            I1 = 0.0;
            I2 = 0.0;

            /**  Gauss - Kronrod 15 - 7 */
            if(order == 15){
                for(i=0;i<15;i++){
                    f_eval = f((b-barc)*KRONROD_NODE_15[i]/2.0+(barc+b)/2.0,param);
                    I2+= f_eval*KRONROD_WEIGHT_15[i];
                    if(i%2 ==1){
                        I1+= f_eval*GAUSS_WEIGHT_7[(i-1)/2];
                    }
                }
            }
            I2 = I2*(b-barc)/2.0;
            I1 = I1*(b-barc)/2.0;
            (*temporaire).I1 = I1;
            (*temporaire).I2 = I2;
        }
        if((*As)->next !=NULL){
            (*As) = (*As)->next;
        }
        else{
            counter = 1;
        }
    };
    return 0;
}

int refine_grid_1D_void(refinment_1D **As, double (*f)(double, void *), void *param, int loga, int order){

    rewind_str_1D(As);
    //printf("After rewind\n");
    int counter = 0;
    while(counter != 1){
        if((*As)->refine == 1){

            //printf("refine_grid_1D 1\n");
            /// Refinment is needed.
            double I1 = 0.0;
            double I2 = 0.0;
            double f_eval;
            int i;

            double a, b, barc;
            a = (*As)->a;
            b = (*As)->b;

            refinment_1D *temporaire = NULL;
            temporaire = (refinment_1D*) malloc(sizeof(refinment_1D));
            //printf("refine_grid_1D 2\n");
            //getchar();

            (*temporaire).refine = 0;
            //printf("refine_grid_1D 2a\n");
            (*temporaire).refineI2g = 1;
            //printf("refine_grid_1D 2b\n");
            (*As)->refineI2g = 1;
             //printf("refine_grid_1D 2c\n");
            (*temporaire).previous = (*As);
            (*temporaire).next = (*As)->next;
            (*As)->next = temporaire;
            //printf("refine_grid_1D 2d\n");
            //getchar();
            if((*temporaire).next !=NULL){
                (*((*temporaire).next)).previous =temporaire;
            }
            //printf("refine_grid_1D 3\n");


            /// proceed to the refinment itself now:
            /// 1) for As:
            if(loga == 0){
                (*As)->a = a;
                barc = (a+b)/2.0;
                (*As)->b = barc;
            }
            else{
                (*As)->a = a;
                barc = exp((log(b)+log(a))/2.0);
                (*As)->b = barc;
            }


            I1 = 0.0;
            I2 = 0.0;
            /**  Gauss - Kronrod 15 - 7 */
            if(order == 15){
                for(i=0;i<15;i++){
                    f_eval = f((barc-a)*KRONROD_NODE_15[i]/2.0+(a+barc)/2.0,param);
                    I2+= f_eval*KRONROD_WEIGHT_15[i];
                    if(i%2 ==1){
                        I1+= f_eval*GAUSS_WEIGHT_7[(i-1)/2];
                    }
                }
            }

            I2 = I2*(barc-a)/2.0;
            I1 = I1*(barc-a)/2.0;
            (*As)->I1 = I1;
            (*As)->I2 = I2;

            /// for temporaire
            (*temporaire).a = barc;            /// The log or non log refinement is already taking care in the c, computed above
            (*temporaire).b = b;

            I1 = 0.0;
            I2 = 0.0;

            /**  Gauss - Kronrod 15 - 7 */
            if(order == 15){
                for(i=0;i<15;i++){
                    f_eval = f((b-barc)*KRONROD_NODE_15[i]/2.0+(barc+b)/2.0,param);
                    I2+= f_eval*KRONROD_WEIGHT_15[i];
                    if(i%2 ==1){
                        I1+= f_eval*GAUSS_WEIGHT_7[(i-1)/2];
                    }
                }
            }
            I2 = I2*(b-barc)/2.0;
            I1 = I1*(b-barc)/2.0;
            (*temporaire).I1 = I1;
            (*temporaire).I2 = I2;
        }
        if((*As)->next !=NULL){
            (*As) = (*As)->next;
        }
        else{
            counter = 1;
        }
    };
    return 0;
}

int free_grid_1D(refinment_1D **As){
    while((*As)->next != NULL){
        (*As) = (*As)->next;
    }
    while((*As)->previous != NULL){
        (*As) = (*As)->previous;
        free((*As)->next);
    }
    return 0;
}


double gauss_kronrod_7_15(double a, double b, double prec, double (*f)(double, double *), double *param, int loga ){


    refinment_1D *refin;
    refin = (refinment_1D *) malloc(1*sizeof(refinment_1D));
    (*refin).a = a;
    (*refin).b = b;
    (*refin).previous = NULL;
    (*refin).next = NULL;


    int i,j;
    double I1 = 0.0;
    double I2 = 0.0;
    double f_eval;


    I1 = 0.0;
    I2 = 0.0;
    for(i=0;i<15;i++){
        f_eval = f((b-a)*KRONROD_NODE_15[i]/2.0+(a+b)/2.0,param);
        I2+= f_eval*KRONROD_WEIGHT_15[i];
        if(i%2 ==1){
            I1+= f_eval*GAUSS_WEIGHT_7[(i-1)/2];
        }
    }
    I2 = I2*(b-a)/2.0;
    I1 = I1*(b-a)/2.0;

    (*refin).I1 = I1;
    (*refin).I2 = I2;
    (*refin).I2g = 0.0;


    j=0;
    int nb_sub;
    //printf("Here err = %g\n", vabs( (I1 - I2 )/I1));

    while ( vabs( (I1 - I2 )/I1) > prec && j<100){
        //printf("Inside 1\n");
        /// Count the elements to get ready for labeling the one which require further subdivision
        nb_sub = Counting_1D( &refin);
        //printf("Inside 2 nb_sub = %d\n", nb_sub);
        rewind_str_1D(&refin);
        //printf("Inside 3\n");
        /// Label the interval which needs refinement
        if(  vabs( ((*refin).I1 - (*refin).I2 )/(*refin).I2) > prec/((double) nb_sub)){
            //printf("Inside 3a\n");
            (*refin).refine = 1;
        }
        //printf("Inside 3aa\n");
        while((*refin).next != NULL){
            //printf("Inside 3b\n");
            refin = (*refin).next;
            if(  vabs( ((*refin).I1 - (*refin).I2 )/(*refin).I2) > prec/((double) nb_sub)){
                (*refin).refine = 1;
            }
        }
        //printf("Inside 4\n");
        /// Proceed with the refinement
        refine_grid_1D(&refin, f, param, loga, 15);
        //printf("Inside 5\n");
        /// Compute the final integrals
        I1 = 0.0;
        I2 = 0.0;
        rewind_str_1D(&refin);
        I1 = (*refin).I1;
        I2 = (*refin).I2;
        while((*refin).next != NULL){
            refin = (*refin).next;
            I1+=(*refin).I1;
            I2+=(*refin).I2;
        }
        //printf("j = %d\n", j);
        //getchar();
        j+=1;
    }
    nb_sub = Counting_1D( &refin);
    //printf("After counting \n");
    //printf("j = %d\t nb_sub= %d\n",j, nb_sub);
    free_grid_1D(&refin);
    //printf("After freeing the grid\n");
    return I2;
}




typedef struct particle{


    /** Energy grid */
    double *E;
    double *Eb;

    /** Distribution function, or coefficients */
    double *f;
    double *f1;
    double *f2;
    double *fderivative;
    double *fderivative1;
    double *fderivative2;

    /** Characteristics of the particle type */
    int dim;
    double emin;
    double emax;

    double mass;
}particle;

int make_particle(particle *pt, double emin, double emax, int dim, double mass){

    int i;
    double step = exp(log(emax/emin)/((double) dim));
    double sqrt_step = sqrt(step);

    pt->emin = emin;
    pt->emax = emax;
    pt->dim = dim;

    /** Energy grid */
    pt->E = calloc(dim, sizeof(double));
    pt->Eb = calloc(dim+1, sizeof(double));

    /** Distribution function, or coefficients */
    pt->f = calloc(dim, sizeof(double));
    pt->f1 = calloc(dim, sizeof(double));
    pt->f2 = calloc(dim, sizeof(double));
    pt->fderivative = calloc(dim, sizeof(double));
    pt->fderivative1 = calloc(dim, sizeof(double));
    pt->fderivative2 = calloc(dim, sizeof(double));

    pt->mass = mass;
    pt->Eb[0] = emin/sqrt_step;
    for(i=0;i<dim;i++){
        pt->E[i] = emin*pow(step,(double) i);
        pt->Eb[i+1] = pt->E[i]*sqrt_step;
    }


    return 0;
};

int remove_grid(particle *pt){

    free(pt->E);
    free(pt->Eb);
    free(pt->f);
    free(pt->fderivative);
    free(pt->fderivative1);
    free(pt->fderivative2);
    return 0;

}


double Int_bpl_dermer_2008(double x, double *param){
    //double K = param[0];
    double p = param[0];
    double gb = param[1];
    if(x < gb){
        return pow(x/gb,-p);
    }
    else{
        return pow(x/gb,(-p-1));
    }
}

/** There is no normalization in that function */
int init_broken_pl_Dermer_2008(particle *pt, double gmin, double gb, double gmax, double p1, double K){


    int i;
    double *param;
    param = calloc(2, sizeof(double));
    //printf("Inside logparabola \tdim = %d\n", pt->dim);
    for(i=0; i<pt->dim; i++){
        //printf(" i = %d\n", i);
        if(pt->Eb[i+1] > gmin && pt->Eb[i]< gmax){
            /** Need integral */
            //printf("LOG PARABOLA  Here 2\n");
            //param[0] = K;
            param[0] = p1;
            param[1] = gb;
            double Imin, Imax;
            Imin =  pt->Eb[i];
            Imax =  pt->Eb[i+1];
            if(gmin > pt->Eb[i] && gmin < pt->Eb[i+1]){
                Imin = gmin;
            }
            if(gmax > pt->Eb[i] && gmax < pt->Eb[i+1]){
                Imax = gmax;
            }
            pt->f[i] =K*gauss_kronrod_7_15(Imin, Imax, 1e-7, &Int_bpl_dermer_2008, param, 0 );

        }
        else{
            pt->f[i]= 0.0;
        }
    }

    double summ = 0.0;
    for(i =0; i<pt->dim;i++){
        summ += pt->f[i];
    }
    for(i =0; i<pt->dim;i++){
        pt->f[i] = pt->f[i]*K/summ;
    }
    free(param);
    return 0.0;
}





/**  ******************
  *  ******************
  *  Synchrotron functions
  *
  */
/** Renormalized emissivities, see notes*/
double Synchrotron_emissivity_photon(double y, double g, double *param){

    double ythreeg2 = y/(3.0*g*g);
    if(ythreeg2 < 150.0){

        double K43 = gsl_sf_bessel_Knu(4.0/3.0, ythreeg2);
        double K13 = gsl_sf_bessel_Knu(1.0/3.0, ythreeg2);
        double first = K43*K13;
        //double second = 0.6*ythreeg2*(K43*K43-K13*K13);
        double second = 0.6*ythreeg2*(K43*K43-K13*K13);
        //printf("ythreeg2 = %g\tresult = %g\t K43 = %g\t K13 = %g\n",ythreeg2, ythreeg2*ythreeg2*(first - second), K43, K13);
        if((first - second) <0.0){
            printf("Function smaller than 0 y = %g\tg = %g\ty/3g2 = %g\tresult = %g\tK43 = %g\t K13 =%g\n", y, g,ythreeg2, (first - second), K43, K13);
            getchar();
        }
        return ythreeg2*ythreeg2*(first - second)/y;
    }
    else{
        return 0;
    }
}

double Synchrotron_emissivity(double y, double g, double *param){

    double ythreeg2 = y/(3.0*g*g);
    if(ythreeg2 < 150.0){

        double K43 = gsl_sf_bessel_Knu(4.0/3.0, ythreeg2);
        double K13 = gsl_sf_bessel_Knu(1.0/3.0, ythreeg2);
        double first = K43*K13;
        //double second = 0.6*ythreeg2*(K43*K43-K13*K13);
        double second = 0.6*ythreeg2*(K43*K43-K13*K13);
        //printf("ythreeg2 = %g\tresult = %g\t K43 = %g\t K13 = %g\n",ythreeg2, ythreeg2*ythreeg2*(first - second), K43, K13);
        if((first - second) <0.0){
            printf("Function smaller than 0 y = %g\tg = %g\ty/3g2 = %g\tresult = %g\tK43 = %g\t K13 =%g\n", y, g,ythreeg2, (first - second), K43, K13);
            getchar();
        }
        return ythreeg2*ythreeg2*(first - second);
    }
    else{
        return 0;
    }
}

/** Expression from Finke, Dermer and Bottcher 2008 */
/** y = nu / nu_B */
double Synchrotron_emissivity_dermer(double y, double g, double *param){

    //double x = (2.0/3.0)*y/(g*g);
    double x = y;
    double logx = log(x)/log(10.0);
    double A0 = x < 1.0 ? -0.35775237 : -0.35842494;
    double A1 = x < 1.0 ? -0.83695385 : -0.79652041;
    double A2 = x < 1.0 ? -1.1449608  : -1.6113032;
    double A3 = x < 1.0 ? -0.68137283 :  0.26055213;
    double A4 = x < 1.0 ? -0.22754737 : -1.6979017;
    double A5 = x < 1.0 ? -0.031967334:  0.032955035;

    double result;
    if(x>=1e-2 && x <= 10){
        result = pow( 10.0, A0 + A1*logx + A2*logx*logx + A3 * pow(logx,3.0) + A4*pow(logx,4.0) + A5*pow(logx,5.0)  );
    }
    else if(x<1e-2){
        result = 1.80842*pow(x,1.0/3.0);
    }
    else{
        result = 0.5*pi*exp(-x)*(1.0-99.0/(162.0*x));
    }
    return result;
}


double Synch_dermer(double x){

    double logx = log(x)/log(10.0);

    double A0 = x < 1.0 ? -0.35775237 : -0.35842494;
    double A1 = x < 1.0 ? -0.83695385 : -0.79652041;
    double A2 = x < 1.0 ? -1.1449608  : -1.6113032;
    double A3 = x < 1.0 ? -0.68137283 :  0.26055213;
    double A4 = x < 1.0 ? -0.22754737 : -1.6979017;
    double A5 = x < 1.0 ? -0.031967334:  0.032955035;

    double result;
    if(x>=1e-2 && x <= 10){
        result = pow( 10.0, A0 + A1*logx + A2*logx*logx + A3 * pow(logx,3.0) + A4*pow(logx,4.0) + A5*pow(logx,5.0)  );
    }
    else if(x<1e-2){
        result = 1.80842*pow(x,1.0/3.0);
    }
    else{
        result = 0.5*pi*exp(-x)*(1.0-99.0/(162.0*x));
    }
    return result;
}


/** They are really the same. */
int comparison_synchrotron_emissivity(){

    FILE *output;
    output = fopen("sync.txt","a");
    double x;
    double step = exp(log(1e15)/100);
    int i;
    double g = 100.0;
    double B= 0.01;
    double nuB = 2.7992483604657531e6*B;
    for(i = 0; i< 100 ; i++){
        x = 1e4*pow(step,i);
        fprintf( output, "%g\t%g\t%g\n", x, 4.0*pi*sqrt3*qe*qe*nuB*Synchrotron_emissivity(x/nuB,  g, &g)/c, sqrt3*pow(qe,3.0)*B*Synchrotron_emissivity_dermer( 4.0*pi*x*c*me/(3.0*qe*B*g*g),  g, &g)/h   );
    }
    fclose(output);
    return 0;

}
/** Function to be integrated
  * *param = x
  */
/*
double Integral_emissivity(double g, double *param){
    return Synchrotron_emissivity(*param , g);
}
*/

/**  3rd of March 2018
  *  Computation of the opacity
  */
double diff_sync_emissivity(double y, double g){

    double x = y/(3.0*g*g);
    if(x < 150.0){

        double K43 = gsl_sf_bessel_Knu(4.0/3.0, x);
        double K13 = gsl_sf_bessel_Knu(1.0/3.0, x);
        double first = K43*K13;
        double second = 1.2*x*(K43*K43-K13*K13);
        //printf("ythreeg2 = %g\tresult = %g\t K43 = %g\t K13 = %g\n",ythreeg2, ythreeg2*ythreeg2*(first - second), K43, K13);
        return -x*y*(2.0/(3.0*g*g*g))*(first - second);
    }
    else{
        return 0.0;
    }

}


double Integral_opacity(double y, double g, double *param){
    /** \\ OLD FORMULATION \\ */
    // double result = (2.0*g/sqrt(g*g-1.0)*Synchrotron_emissivity( y, g, &y) + diff_sync_emissivity( y, g))/(y*y);
    /** \\ NEW FORMULATION \\ */
    double result = ((2.0*g*g-1.0)/(g*(g*g-1.0))*Synchrotron_emissivity( y, g, &y) + diff_sync_emissivity( y, g))/(y*y);
    return result;
    //Synchrotron_emissivity(*param , g, param);
}




/**  3rd of March 2018
  *  Computation of the emissivity
  */
int sync_emissivity(double **syncrate, particle *ph, particle *elec){

    /**  If the file exist, just read.
      */

    int i,j;
    double a;
    FILE *intput;
    intput = fopen("synchrotron_rate.txt", "r");
    while(fscanf(intput, "%d%d%lf", &i,&j,&a) != EOF){
        syncrate[i][j] = a;
    }
    fclose(intput);
    printf("Sync emissivity loaded \n");
    //getchar();

    /**  If not, create */
    /*
    int i, j;
    double *a, *b;
    //double param;
    a = (double *) calloc(2, sizeof(double));
    b = (double *) calloc(2, sizeof(double));

    FILE *output;
    output = fopen("synchrotron_rate.txt", "a");
    //output1 = fopen("synchrotron_rate_gnuplot.txt", "a");

    double result;
    for(i=0;i<GRID_P;i++){
        a[0] = ph->Eb[i];
        a[1] = ph->Eb[i+1];;

        for(j=0;j<GRID_E;j++){
            b[0] = elec->Eb[j];
            b[1] = elec->Eb[j+1];
            //printf("i = %d\t j = %d\tresult = %g\tx = %g\t%g\tg = %g\t%g\n",i,j, result, a[0],a[1], b[0], b[1]);
            result = gauss_kronrod_7_15_multi2_gauss(a, b, 1e-4, &Synchrotron_emissivity_photon, a, 0, "integral.txt");
            //result = gauss_kronrod_7_15( a[0], a[1], 1e-4, &Integral_emissivity, &param, 0 );

            syncrate[i][j] = result;

            //printf(" i = %d\t j = %d\t x = %g - %g\t g = %g - %g\tresult = %g\n", i, j, a[0], a[1], b[0], b[1], result);
            fprintf(output, "%d\t%d\t%g\n",i,j,result);
            //fprintf(output1, "%d\t%d\t%g\n",i,j,result);
            //getchar();
        }
        //fprintf(output1, "\n");

        //getchar();
    }
    fclose(output);
    free(a);
    free(b);
    */
    return 0;

}



/**  3rd of March 2018
  *  Computation of the opacity
  */
int opacity(double **opacity_rate, particle *ph, particle *elec){

    /**  If the file exist, just read.
      */

    int i,j;
    double a;
    FILE *intput;
    intput = fopen("tau.txt", "r");
    while(fscanf(intput, "%d%d%lf", &i,&j,&a) != EOF){
        opacity_rate[i][j] = a;
    }
    fclose(intput);

    /**  If not, create */
    /*
    int i, j;
    double *a, *b;
    //double param;
    a = (double *) calloc(2, sizeof(double));
    b = (double *) calloc(2, sizeof(double));

    FILE *output;
    //FILE *output1;
    output = fopen("tau.txt", "a");
    // output1 = fopen("tau_gnuplot.txt", "a");

    double result;
    for(i=0;i<GRID_P;i++){
        a[0] = ph->Eb[i];
        a[1] = ph->Eb[i+1];
        for(j=0;j<GRID_E;j++){
            b[0] = elec->Eb[j];
            b[1] = elec->Eb[j+1];
            //param = ph->E[i];
            printf("i = %d\t j = %d\tresult = %g\n",i,j, result);
            result = gauss_kronrod_7_15_multi2_gauss(a, b, 1e-4, &Integral_opacity, a, 0, "integral.txt");
            //result = gauss_kronrod_7_15( a[0], a[1], 1e-4, &Integral_opacity, &param, 0 );

            opacity_rate[i][j] = result;
            //printf(" i = %d\t j = %d\t x = %g - %g\t g = %g - %g\tresult = %g\n", i, j, a[0], a[1], b[0], b[1], result);
            fprintf(output, "%d\t%d\t%g\n",i,j,result);
      //      fprintf(output1, "%d\t%d\t%g\n",i,j,result);
            //getchar();
        }
        //fprintf(output1, "\n");
        //getchar();
    }
    fclose(output);
    //fclose(output1);
    //printf("result = %g\n",result);
    free(a);
    free(b);
    */
    return 0;
}


/**  ************************************ */
/**  ************************************ */
/**  ************************************ */
/**  4/03/2018
  *  Pair kernel.
  */
double phipair(double x){
    /// OLD
    /*
    double v = x -1.0;
    double sqrtv = sqrt(v);
    double sqrtx = sqrt(x);
    double w = x < 1e8 ? (sqrtx+sqrtv)/(sqrtx-sqrtv) : 4.0*x;
    double wp1 = w+1.0;
    double logw = log(w);
    double logw2 = logw*logw;
    double logwp1 = log(w+1.0);
    double logwp12 = logwp1*logwp1;
    double dilog = gsl_sf_dilog(1.0/wp1);
    //printf("x = %g\t dilog = %g\n",x, dilog);
    //if( (v!=v)  || (w!=w) || (logw != logw) || (logwp12 != logwp12) || (dilog != dilog))
    //    printf("nan:  x = %g\t v = %g\t w = %g\t logw = %g\t logw12 = %g\t dilog = %g\n",x, v, w, logw, logwp12, dilog);
    //    getchar();


    //
    //printf("nan:  x = %g\t v = %g\t w = %g\t logw = %g\t logw12 = %g\t dilog = %g\n",x, v, w, logw, logwp12, dilog);
    return (2.0*v*x+1.0)*logw/x - 2.0*(v+x)*sqrtv/sqrtx - logw2 + 2.0*logwp12 + 4.0 * dilog  - pi2/3.0;
    */
    /// NEW
    double v = x -1.0;
    double w = x < 1e8 ? (sqrt(v+1.0) + sqrt(v))/(sqrt(v+1)-sqrt(v)) : 4.0*x;
    double sqrtv = sqrt(v);
    double sqrtvp1 = sqrt(v+1.0);
    double wp1 = w+1.0;
    double logw = log(w);
    double logw2 = logw*logw;
    double logwp1 = log(w+1.0);
    double logwp12 = logwp1*logwp1;
    double dilog = gsl_sf_dilog(1.0/wp1);
    printf("dilog = %g\n",dilog);
    //printf("x = %g\t dilog = %g\n",x, dilog);
    //if( (v!=v)  || (w!=w) || (logw != logw) || (logwp12 != logwp12) || (dilog != dilog))
    //    printf("nan:  x = %g\t v = %g\t w = %g\t logw = %g\t logw12 = %g\t dilog = %g\n",x, v, w, logw, logwp12, dilog);
    //    getchar();


    //
    //printf("nan:  x = %g\t v = %g\t w = %g\t logw = %g\t logw12 = %g\t dilog = %g\n",x, v, w, logw, logwp12, dilog);
    // printf("log(10 = %g\n",log(10));
    return (2.0*v*v+2.0*v+1.0)*logw/(v+1.0)  - 2.0*(2.0*v+1.0)*sqrtv/sqrtvp1 - logw2 + 2.0*logwp12 + 4.0*dilog - pi2/3.0;

}

double Rpp(double x, double *param){
    if(x>1.0){
        return phipair(x)/(x*x);
    }
    else{
        return 0.0;
    }
}

int test_Rpp(){
    double a = 2.0;
    printf("Rpp (50) = %g\n", Rpp(50,&a));
    printf("Rpp (10) = %g\n", Rpp(10,&a));
    printf("Rpp (5) = %g\n", Rpp(5,&a));
    return 0;
}

int pair_opacity(double **opacity_rate, particle *ph){

    /**  If the file exist, just read.
      */

    int i,k;
    double a;
    FILE *intput;
    intput = fopen("pair_rate.txt", "r");
    while(fscanf(intput, "%d%d%lf", &i,&k,&a) != EOF){
        opacity_rate[i][k] = a;
    }
    fclose(intput);


    /**  If not, create */
    /*

    int i, k;


    FILE *output;
    output = fopen("pair_rate.txt", "a");

    double result;
    for(i=0;i<GRID_P;i++){                /// outgoing photon
        for(k = 0;k<GRID_P;k++){            /// incoming photon

   //i = 101;
   //k = 107;

            result = gauss_kronrod_7_15(  (ph->Eb[k]*h/mec2)*(ph->E[i]*h/mec2), (ph->Eb[k+1]*h/mec2)*(ph->E[i]*h/mec2), 1e-4, &Rpp, &result, 0 );
            opacity_rate[i][k] = (3.0*sigmaT)*result/8.0;
   //printf( "low = %g\t up = %g\t result = %g\t ", (ph->Eb[k]*h/mec2)*(ph->E[i]*h/mec2), (ph->Eb[k+1]*h/mec2)*(ph->E[i]*h/mec2), (3.0*sigmaT)*result/8.0 );
   //getchar();
            fprintf(output, "%d\t%d\t%g\n", i, k, (3.0*sigmaT)*result/8.0);
        }
    }

    fclose(output);
    */
    return 0;
}


double Jones_head_on(double x1, double x, double g){   /// x_1 to x
    /*
    Jones provided an asymptotic expression for compton.
    */
    double q = x/(4.0*x1*g*g*(1.0-x/g));
    if(q>0.25/g/g && q<=1.0){
        double fqG = 2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+0.5*(16.0*x1*g*q*x1*g*q)*(1.0-q)/(1.0+4.0*x1*g*q);
        return fqG/(x1*g*g);
    }
    else{
        return 0.0;
    }

}

double Integral_Compton(double *x, size_t dim, void *params){
    // return x[0]*Jones_head_on(x[1], x[0], x[2]);
    return Jones_head_on(x[1], x[0], x[2]);
}

int inverse_compton(double ***opacity_rate, particle *ph, particle *elec){

    /**  If the file exist, just read.
      */

    int i,j,k;
    double a,err;
    FILE *intput;
    intput = fopen("IC_rate.txt", "r");
    while(fscanf(intput, "%d%d%d%lf%lf", &i,&j,&k,&a,&err) != EOF){
        opacity_rate[i][j][k] = a;
    }
    fclose(intput);


    /**  If not, create */
    /*
    int i, j, k;


    FILE *output;
    //FILE *output1;

    double xl[3];
    double xu[3];

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &Integral_Compton, 3, 0 };

    size_t calls = 50000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);




    output = fopen("IC_rate.txt", "a");
    // output1 = fopen("IC_rate_gnuplot.txt", "a");


    double result, err;
    for(i=0;i<GRID_P;i++){                /// outgoing photon
        xl[0] = ph->Eb[i]*h/mec2;
        xu[0] = ph->Eb[i+1]*h/mec2;
        for(k = 0;k<i+1;k++){            /// incoming photon
            printf("i = %d\t k =%d\n", i, k);
            xl[1] = ph->Eb[k]*h/mec2;
            xu[1] = ph->Eb[k+1]*h/mec2;
            for(j=0;j<GRID_E;j++){
                xl[2] = elec->Eb[j];
                xu[2] = elec->Eb[j+1];
                //result = gauss_kronrod_7_15_multi2_gauss(a, b, 1e-4, &Integral_opacity, a, 0, "integral.txt");
                gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s, &result, &err);
                //result = gauss_kronrod_7_15( a[0], a[1], 1e-4, &Integral_opacity, &param, 0 );

                opacity_rate[i][k][j] = result;
                // printf("i = %d\t Ei = %g\t k = %d\t Ek = %g\t  j = %d\t Ej = %g\tresult = %g\n", i, ph->E[i]*h/mec2, k,ph->E[k]*h/mec2,j, elec->E[j], result);
                //printf(" i = %d\t j = %d\t x = %g - %g\t g = %g - %g\tresult = %g\n", i, j, a[0], a[1], b[0], b[1], result);
                fprintf(output, "%d\t%d\t%d\t%g\t%g\n",i,k,j,result,err);
                //fprintf(output1, "%d\t%d\t%d\t%g\n",i,k,j,result);
                //getchar();
            }
            //fprintf(output1, "\n");
        //getchar();
        }
    }
    gsl_monte_plain_free (s);
    fclose(output);
    //fclose(output1);
    printf("result = %g\n",result);
    */
    return 0;
}


int rescal_grid_old(particle ph1, particle ph2, double nuB){
    /** Rescal the energy from 1 grid to the other */
    double nubar = ph1.Eb[0]*nuB;
    int i,j;
    /** Case A : nubar < ph2.Eb[0], that is to say we have to forget some of the points*/
    if(nubar < ph2.Eb[0]){
        i = 0;
        printf("Case A \n");
        while(nubar < ph2.Eb[0]){
            i += 1;
            nubar = ph1.Eb[i]*nuB;
        }
        i = i-1; /// Now  ph1.Eb[i] < ph2.Eb[0] <= ph1.Eb[i+1]
        for(j = i;j< GRID_P; j++){
            ph2.f[j-i] = ph1.f[j]*( ph1.Eb[j+1]*nuB - ph2.Eb[j-i])/((ph1.Eb[j+1] - ph1.Eb[j])*nuB) + ph1.f[j+1]*( ph2.Eb[j-i+1] - ph1.Eb[j+1]*nuB )/((ph1.Eb[j+2] - ph1.Eb[j+1])*nuB);
        }
    }
    else{       /** Case B : nubar > ph2.Eb[0] */
        i = 0;
        printf("Case B \n");
        while(nubar > ph2.Eb[i]){
            i += 1;
        }
        i = i-1; /// Now ph2.Eb[i] < ph1.Eb[0] * nuB < ph2.Eb[i+1]
        j = i;
        ph2.f[j] = ph1.f[j-i]*( ph2.Eb[j+1] - ph1.Eb[j-i]*nuB )/((ph1.Eb[j-i+1] - ph1.Eb[j-i])*nuB);
        for(j = i+1;j< GRID_P-1; j++){
            ph2.f[j] = ph1.f[j-i]*( ph2.Eb[j+1] - ph1.Eb[j-i]*nuB )/((ph1.Eb[j-i+1] - ph1.Eb[j-i])*nuB) + ph1.f[j-i-1]*(- ph2.Eb[j] + ph1.Eb[j-i]*nuB )/((ph1.Eb[j-i] - ph1.Eb[j-i-1])*nuB);
        }
    }
    return 0;
}


int rescal_grid(particle ph1, particle ph2, double nuB){
    /** Rescal the energy from 1 grid to the other */
    double nubar = ph1.Eb[0]*nuB;
    int i,j;
    /** Case A : nubar < ph2.Eb[0], that is to say we have to forget some of the points*/
    if(nubar < ph2.Eb[0]){
        i = 0;
        //printf("Case A \n");
        while(nubar < ph2.Eb[0]){
            i += 1;
            nubar = ph1.Eb[i]*nuB;
        }
        i = i-1; /// Now  ph1.Eb[i] < ph2.Eb[0] <= ph1.Eb[i+1]
        for(j = i;j< GRID_P; j++){
            ph2.f[j-i] = (ph1.f[j]*( ph1.Eb[j+1]*nuB - ph2.Eb[j-i]) + ph1.f[j+1]*( ph2.Eb[j-i+1] - ph1.Eb[j+1]*nuB ))/(ph2.Eb[j-i+1] - ph2.Eb[j-i]);

        }
    }
    else{       /** Case B : nubar > ph2.Eb[0] */
        i = 0;
        printf("Case B \n");
        while(nubar > ph2.Eb[i]){
            i += 1;
        }
        i = i-1; /// Now ph2.Eb[i] < ph1.Eb[0] * nuB < ph2.Eb[i+1]
        j = i;
        ph2.f[j] = ph1.f[j-i]*( ph2.Eb[j+1] - ph1.Eb[j-i]*nuB )/(ph2.Eb[j+1] - ph2.Eb[j]);

        for(j = i+1;j< GRID_P-1; j++){
            ph2.f[j] = (ph1.f[j-i]*( ph2.Eb[j+1] - ph1.Eb[j-i]*nuB ) + ph1.f[j-i-1]*(- ph2.Eb[j] + ph1.Eb[j-i]*nuB ))/(ph2.Eb[j+1] - ph2.Eb[j]);
        }
    }
    return 0;
}


double tau_EBL(double nu ){

    double A = 0.010677;
    double B = -0.238895;
    double C = -1.004;
    double D = 54.1465;
    double E = -313.486;

    double EeV = log(nu*Hz_to_eV)/log(10.0);
    double tau = exp(A*pow(EeV,4.0)+ B*pow(EeV,3.0)+C*pow(EeV,2.0)+D*EeV + E);
    //printf("tau = %g\t EeV = %g\targ = %g\n", tau, EeV,A*pow(EeV,4.0)+ B*pow(EeV,3.0)+C*pow(EeV,2.0)+D*EeV + E);
    double result = 1.0;
    if(EeV > 10){
        result = (1.0-exp(-tau))/tau;
    }
    return result;

}

int Emissivity(double Bfield, double delta, double R, double d, particle photon_IC, particle electron){

    /** Bfield is the comoving magnetic field */
    double delta3 = delta*delta*delta;

    double nuB = 2.7992483604657531e6*Bfield;
    double A = 4.0*pi*sqrt3*qe*qe/c/h;

    double ASSA = 4.0*pi*sqrt3*qe*qe*nuB/c;

    double tauCst = R/(8.0*pi*me*nuB*nuB);


    particle photon;
    double emin = 0.01;
    double emax = 1e30;
    make_particle(&photon, emin, emax, GRID_P, 0.0);        /// Units of nu_B

    /** Load the cross-section and rates */
    int i,j,k;
    double **syncrate;
    syncrate = (double **) malloc(GRID_P*sizeof(double *));
    for(i = 0 ; i< GRID_P; i++){
        syncrate[i] = (double *) calloc(GRID_E, sizeof(double));
    }
    sync_emissivity(syncrate, &photon, &electron);


    double **opacity_rate;
    opacity_rate = (double **) malloc(GRID_P*sizeof(double *));
    for(i = 0 ; i< GRID_P; i++){
        opacity_rate[i] = (double *) calloc(GRID_E, sizeof(double));
    }
    opacity(opacity_rate, &photon, &electron);


    double ***IC_rate;
    IC_rate = (double ***) malloc(GRID_P*sizeof(double **));
    for(i = 0 ; i< GRID_P; i++){
        IC_rate[i] = (double **) malloc(GRID_P* sizeof(double*));
        for(j = 0; j<GRID_P; j++){
            IC_rate[i][j] = (double *) calloc(GRID_E, sizeof(double));
        }
    }
    inverse_compton(IC_rate, &photon, &electron);
/*
    double **pair_rate;
    pair_rate = (double **) malloc(GRID_P*sizeof(double *));
    for(i = 0 ; i< GRID_P; i++){
        pair_rate[i] = (double *) calloc(GRID_P, sizeof(double));
    }
*/


    /** Begining of the computation */
    clock_t begin = clock();

    /** Synchrotron */
    /** 15/05/2018: this is working, comparison was made with Tramacere et al. on their web application*/
    for(i = 0; i< GRID_P;i++){
        double summ = 0.0;
        for(j = 0;j< GRID_E;j++){
            summ += A*syncrate[i][j]*(electron.f[j]/(electron.Eb[j+1]-electron.Eb[j]));
        }
        photon.f[i] = summ/(photon.Eb[i+1]-photon.Eb[i]);
        //printf(" i = %d\t photon i = %g\n", i , photon.f[i]);
    }

    /**  SSA */

    for(i = 0; i <GRID_P ;i++){
        double summtau = 0.0;
        for(j = 0;j< GRID_E;j++){
            summtau += opacity_rate[i][j]*(electron.f[j]/(electron.Eb[j+1]-electron.Eb[j]));
        }

        summtau = ASSA*tauCst*summtau;

        //photon.f[i] = photon.f[i]*onemexpdexp(summtau);
    }

    rescal_grid(photon, photon_IC, nuB);

    /** Count photons */
    /*double summB = 0.0;
    double summA = 0.0;
    for(i = 0; i < GRID_P; i++){
        summB += photon.f[i] * (photon.Eb[i+1]-photon.Eb[i])*nuB;
        summA += photon_IC.f[i] * (photon_IC.Eb[i+1]-photon_IC.Eb[i]);
    }
    printf("Photon number before = %g\t after = %g\n", summB, summA);
    */

    /** IC */
    double summIC = 0.0;
    double summICfinal;
    double summtemp;
    // double summ_ph_IC = 0.0;
    for(i = 0; i <GRID_P ;i++){
        summIC = 0.0;
        for(j = 0;j< i-1;j++){
            summtemp = 0.0;
            for(k = 0; k< GRID_E;k++){
                summtemp += IC_rate[i][j][k]*electron.f[k]/(electron.Eb[k+1]-electron.Eb[k]);
            }
            summIC += photon_IC.f[j]*summtemp;
        }
        summICfinal = 0.75*c*sigmaT*summIC /(photon_IC.Eb[i+1]-photon_IC.Eb[i]); ///(photon_IC.Eb[i+1]-photon_IC.Eb[i]);

        /// Dermer : multiply the emissivity by the volume.
        photon_IC.f2[i] = (4.0*pi*R*R*R*delta3/(d*Mpc_to_cm*d*Mpc_to_cm*4.0*pi*3))*( (R/c)*summICfinal)*h*photon_IC.E[i] ;

    }
    for(i = 0; i < GRID_P; i++){
        photon_IC.f[i] = (4.0*pi*R*R*R*delta3/(d*Mpc_to_cm*d*Mpc_to_cm*4.0*pi*3))*(photon_IC.f[i] )*h*photon_IC.E[i] ;
    }
    /** Pair */
    double summpair;
    double tau_pair;
    for(i = 0; i <GRID_P-1 ;i++){
        summpair = 0.0;
        for(j = 0;j< GRID_P;j++){
        //    summpair += pair_rate[i][j]*photon_IC.f2[j]/(h*photon_IC.E[j]);
        }
        tau_pair = 4.0*pi*summpair*R/(delta*h);
    }


    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent = %g\n", time_spent);

    for(i = 0 ; i< GRID_P; i++){
        free(syncrate[i]);
        //free(opacity_rate[i]);
    }
    free(syncrate);
    //free(opacity_rate);

    remove_grid(&photon);

/*
    for(i = 0 ; i< GRID_P; i++){
        for(j = 0; j<GRID_P; j++){
            free(IC_rate[i][j]);
        }
        free(IC_rate[i]);
    }
    free(IC_rate);
*/

	return 0;
}


/** We use formula 21 of Finke, Dermer and Bottchner to compute the synchrotron spectrum */
int Dermer_exact_paper_computation(){


    int i,j;


    /** Param 3a : D07 IBL tv = 30s*/
    double Bfield = 0.088;
    double delta = 230;
    double R = 0.19e15;
    double gc = 3.1e4;
    double gmax = 1.3e5;
    double Ke = 2e40;
    double V = (4.0*pi*R*R*R/3.0);


    double d = 540.0;   /// Distance in Mpc
    double z = 0.116;

    /** Initialization */
    double emin = 1e-15;
    double emax = 1e15;
    particle photon;
    particle photon_Dermer;
    particle electron;

    make_particle(&photon_Dermer, emin, emax, GRID_P, 0.0);       /// Units of nu
    make_particle(&photon, 0.01, 1e30, GRID_P, 0.0);       /// Units of nu
    make_particle(&electron, 10.0, 1e8, GRID_E,0.0);

    /** Init electrons */
    /// Param 1
    //init_broken_pl_Dermer_2008(&electron, 1e3, 1.1e5, 1.3e5, 2.7, 9e38/(4.0*pi*pow(R,3.0)/3.0));
    /// Param 2
    //init_broken_pl_Dermer_2008(&electron, 1e3, 2.1e5, 9.4e5, 2.7, 3e41/(4.0*pi*pow(R,3.0)/3.0));

    /// Param 3a
    init_broken_pl_Dermer_2008(&electron, 1e3, 3.1e4, 1.3e5, 2.7, 2e40*delta*delta*delta/V);


    double summe = 0.0;
    for(i = 0; i<GRID_E;i++){
        summe += electron.f[i];
    }


    double *param;
    param = (double *) calloc(4, sizeof(double));
    param[1] = Bfield;
    param[2] = gc;
    param[3] = Ke/V;


    Emissivity(Bfield, delta, R,  d, photon, electron);




    FILE *distrib;
    distrib = fopen("Dermer_paper_a.txt","a");


    /** Export data */
    for(i = 0; i <GRID_P ;i++){
        fprintf(distrib,"%g\t%g\t%g\t%g\t%g\n", photon.E[i]*delta/(1.0+z), photon.f[i], photon.f2[i], photon_Dermer.E[i]*delta*(mec2/h)/(1.0+z),photon_Dermer.f2[i] );
    }


    remove_grid(&photon_Dermer);
    remove_grid(&photon);
    remove_grid(&electron);


    return 0;

}

int main(){

    Dermer_exact_paper_computation();

    return 0;
}
