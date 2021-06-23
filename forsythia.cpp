#include "forsythia.hpp"

#include "fp2_arithm.hpp"
#include "ec_arithm.hpp"

ForsythiaParamSet g_param_set[2] = 
{
    // forsythia80
    {
        .p = "2085906331811366794969677127350340075070763705125644383096310746405164434883346431",
        .curve_coeffs =
        {
            .A =
            {
                "1697515660496260333664277950769004218366836292978277269037599020171750501119102405",
                "0"
            },
            .B =
            {
                "1", "0"
            },
            .C =
            {
                "1", "0"
            }
        },
        .torsion_params =
        {
            // for_side_2
            {
                .e = "137",
                .xP =
                {
                    "522422669065078842121178374075221380702967188678242626309221604145431102848928362",
                    "247405703781105551606345687695260444799049769254304627035480764604919092477393513"
                },
                .xQ =
                {
                    "1935719375068975114132323850650997729213850194744870562688286399344274453352521538",
                    "1392414880196036042209526661794485734467392457965531424718263704605987281571078758"
                },
                .xR =
                {
                    "326730451538212787712878139069670336778338106785713016573865175340742213327321",
                    "1537197473986856443425048184105212707493203641760886610280536355975232918588188578"
                }
            },
            // for_side_3
            {
                .e = "84",
                .xP =
                {
                    "1242019673579151378282358573982204479435854333343353419255048183921055246579249581",
                    "1666773317262840167533615209246360419786884541146535296789773358382685934474358259"
                },
                .xQ =
                {
                    "1870312412146385036113006824318063366443144427846925504122731796477963524345584993",
                    "101381397618719655616512673486178812026799883731587251393744313522319210756895963"
                },
                .xR =
                {
                    "1355940430959122699383977573426359021645384327366856460395952015098679236349062962",
                    "2051183253123910664141141590855659995581465394826565067955000590109229293369149536"
                }
            }
        }
    },
    // forsythia128
    {
        .p =
        "72753009203726392079824034218111442195971520406335723496697481258231706650008027805076328969921155924624196919691420358410239",
        .curve_coeffs =
        {
            .A =
            {
                "17493471122392659065800712772515557709670529132138038021277127480141235217844132013314393031224679806017205814419946118575015",
                "0"
            },
            .B =
            {
                "1", "0"
            },
            .C =
            {
                "1", "0"
            }
        },
        .torsion_params =
        {
            // for_side_2
            {
                .e = "208",
                .xP =
                {
                    "22844706363854284894147775134188266748171037718525606388952092973941461889797810760893250673350348263076107008036776087419684",
                    "63841253474393277610681807551511151174539169124225858012406716798726146657651778558584054371003621807894257918476250082064384"
                },
                .xQ =
                {
                    "34044695409104666325981838619655510021087561353322719531337987530201536824741279043514521861509851315864518032563593765265552",
                    "7452959005574056344724828511102403135482654493550936362957164066447002781721529155722396704782067294511069235423691009394527"
                },
                .xR =
                {
                    "25778018728727103180300005328978349553036935310641789949646362225128100781829983091015333508539827341144722183109249620913758",
                    "11413233635764865559824758198216090792700224793521820794264309602886357645919539783830031513188654959267247865901357493441175"
                }
            },
            // for_side_3
            {
                .e = "129",
                .xP =
                {
                    "25887208518686915305880592111505623667305562112885313561494547900787514982371722513125288743654421705575186856690481973538601",
                    "63171582405850212153392615420754749458906109385443911908471795928875505646527631230543450128882144819167248880126972102946985"
                },
                .xQ =
                {
                    "51246862455156643703989113914829868661496743069715159976618019431407888217702674106039330579553028268602859658895210237203430",
                    "65722901043151396079631139580312716370936696398017584560491049904137092648979706025324108129233876087941727320577731828052182"
                },
                .xR =
                {
                    "22045120116694800531460819349229131207455760935946806579883662380096575369678802160807315412792756473501707424823553744928011",
                    "54405859870116959587573379590518177683471769362657147132415860485392441225213187894762521827615463033582873359637289624803801"
                }
            }
        }
    },
};