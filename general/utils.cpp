// ============================================================================
/// \file utils.cpp
// ============================================================================
/// \author Jason Sanders
/// \date 2014-2015
/// Institute of Astronomy, University of Cambridge (and University of Oxford)
// ============================================================================

// ============================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ============================================================================
/// \brief Utility functions
//============================================================================

/*======================================*/
/*        A few utility functions       */
/*======================================*/

#include "utils.h"

const double weights4[4] = {0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763};
const double abscissa4[4] = {0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363};
const double weights5[5] = {0.2955242247147529,0.2692667193099963,0.2190863625159820,0.1494513491505806,0.0666713443086881};
const double abscissa5[5] = {0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717};
const double weights8[8]={0.1894506104550685,0.1826034150449236,0.1691565193950025,0.1495959888165767,0.1246289712555339,0.0951585116824928,0.0622535239386479,0.0271524594117541};
const double abscissa8[8]={0.0950125098376374,0.2816035507792589,0.4580167776572274,0.6178762444026438,0.7554044083550030,0.8656312023878318,0.9445750230732326,0.9894009349916499};
const double abscissa10[10] = {0.0765265211334973337546404,0.2277858511416450780804962,0.3737060887154195606725482,0.5108670019508270980043641,0.6360536807265150254528367,0.7463319064601507926143051,0.8391169718222188233945291,0.9122344282513259058677524,0.9639719272779137912676661,0.9931285991850949247861224};
const double weights10[10] = {0.1527533871307258506980843,0.1491729864726037467878287,0.1420961093183820513292983,0.1316886384491766268984945,0.1181945319615184173123774,0.1019301198172404350367501,0.0832767415767047487247581,0.0626720483341090635695065,0.0406014298003869413310400,0.0176140071391521183118620};
const double abscissa16[16] = {0.0483076656877383,0.1444719615827965,0.2392873622521371,0.3318686022821277,0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816};
const double weights16[16]={0.0965400885147278,0.0956387200792749,0.0938443990808046,.0911738786957639,0.0876520930044038,0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,0.0253920653092621,0.0162743947309057,0.0070186100094701};
const double abscissa32[32]={0.0243502926634244325089558,0.0729931217877990394495429,0.1214628192961205544703765,0.1696444204239928180373136,0.2174236437400070841496487,0.2646871622087674163739642,0.3113228719902109561575127,0.3572201583376681159504426,0.4022701579639916036957668,0.4463660172534640879849477,0.4894031457070529574785263,0.5312794640198945456580139,0.5718956462026340342838781,0.6111553551723932502488530,0.6489654712546573398577612,0.6852363130542332425635584,0.7198818501716108268489402,0.7528199072605318966118638,0.7839723589433414076102205,0.8132653151227975597419233,0.8406292962525803627516915,0.8659993981540928197607834,0.8893154459951141058534040,0.9105221370785028057563807,0.9295691721319395758214902,0.9464113748584028160624815,0.9610087996520537189186141,0.9733268277899109637418535,0.9833362538846259569312993,0.9910133714767443207393824,0.9963401167719552793469245,0.9993050417357721394569056};
const double weights32[32]={0.0486909570091397203833654,0.0485754674415034269347991,0.0483447622348029571697695,0.0479993885964583077281262,0.0475401657148303086622822,0.0469681828162100173253263,0.0462847965813144172959532,0.0454916279274181444797710,0.0445905581637565630601347,0.0435837245293234533768279,0.0424735151236535890073398,0.0412625632426235286101563,0.0399537411327203413866569,0.0385501531786156291289625,0.0370551285402400460404151,0.0354722132568823838106931,0.0338051618371416093915655,0.0320579283548515535854675,0.0302346570724024788679741,0.0283396726142594832275113,0.0263774697150546586716918,0.0243527025687108733381776,0.0222701738083832541592983,0.0201348231535302093723403,0.0179517157756973430850453,0.0157260304760247193219660,0.0134630478967186425980608,0.0111681394601311288185905,0.0088467598263639477230309,0.0065044579689783628561174,0.0041470332605624676352875,0.0017832807216964329472961};

const double abscissa50[50]={0.0156289844215430828722167,0.0468716824215916316149239,0.0780685828134366366948174,0.1091892035800611150034260,0.1402031372361139732075146,0.1710800805386032748875324,0.2017898640957359972360489,0.2323024818449739696495100,0.2625881203715034791689293,0.2926171880384719647375559,0.3223603439005291517224766,0.3517885263724217209723438,0.3808729816246299567633625,0.4095852916783015425288684,0.4378974021720315131089780,0.4657816497733580422492166,0.4932107892081909335693088,0.5201580198817630566468157,0.5465970120650941674679943,0.5725019326213811913168704,0.5978474702471787212648065,0.6226088602037077716041908,0.6467619085141292798326303,0.6702830156031410158025870,0.6931491993558019659486479,0.7153381175730564464599671,0.7368280898020207055124277,0.7575981185197071760356680,0.7776279096494954756275514,0.7968978923903144763895729,0.8153892383391762543939888,0.8330838798884008235429158,0.8499645278795912842933626,0.8660146884971646234107400,0.8812186793850184155733168,0.8955616449707269866985210,0.9090295709825296904671263,0.9216092981453339526669513,0.9332885350430795459243337,0.9440558701362559779627747,0.9539007829254917428493369,0.9628136542558155272936593,0.9707857757637063319308979,0.9778093584869182885537811,0.9838775407060570154961002,0.9889843952429917480044187,0.9931249370374434596520099,0.9962951347331251491861317,0.9984919506395958184001634,0.9997137267734412336782285};
const double weights50[50]={0.0312554234538633569476425,0.0312248842548493577323765,0.0311638356962099067838183,0.0310723374275665165878102,0.0309504788504909882340635,0.0307983790311525904277139,0.0306161865839804484964594,0.0304040795264548200165079,0.0301622651051691449190687,0.0298909795933328309168368,0.0295904880599126425117545,0.0292610841106382766201190,0.0289030896011252031348762,0.0285168543223950979909368,0.0281027556591011733176483,0.0276611982207923882942042,0.0271926134465768801364916,0.0266974591835709626603847,0.0261762192395456763423087,0.0256294029102081160756420,0.0250575444815795897037642,0.0244612027079570527199750,0.0238409602659682059625604,0.0231974231852541216224889,0.0225312202563362727017970,0.0218430024162473863139537,0.0211334421125276415426723,0.0204032326462094327668389,0.0196530874944353058653815,0.0188837396133749045529412,0.0180959407221281166643908,0.0172904605683235824393442,0.0164680861761452126431050,0.0156296210775460027239369,0.0147758845274413017688800,0.0139077107037187726879541,0.0130259478929715422855586,0.0121314576629794974077448,0.0112251140231859771172216,0.0103078025748689695857821,0.0093804196536944579514182,0.0084438714696689714026208,0.0074990732554647115788287,0.0065469484508453227641521,0.0055884280038655151572119,0.0046244500634221193510958,0.0036559612013263751823425,0.0026839253715534824194396,0.0017093926535181052395294,0.0007346344905056717304063};