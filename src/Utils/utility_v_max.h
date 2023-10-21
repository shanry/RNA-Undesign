/*
 *utility_v_max.h*
 provides max and min score for vienna energy function.

 author: Wei Yu Tang (based on Kai Zhao, Dezhong Deng code)
 edited by: 10/2023
*/

#ifndef FASTCKY_UTILITY_V_MAX_H
#define FASTCKY_UTILITY_V_MAX_H

// pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
// nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U

#include "utility_v.h"


inline int v_score_hairpin_max(int i, int j, int tetra_hex_tri_index = -1) {
    int nuci = 3, nuci1 = 1, nucj_1 = 1, nucj = 4;
    int size = j-i-1;
    int type = NUC_TO_PAIR(nuci, nucj);

    int energy;

    if(size <= 30)
        energy = hairpin37[size];
    else
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));

    if(size < 3) return energy; /* should only be the case when folding alignments */
// #ifdef SPECIAL_HP
//     // if(special_hp){
//         if (size == 4 && tetra_hex_tri_index > -1)
//             return Tetraloop37[tetra_hex_tri_index];
//         else if (size == 6 && tetra_hex_tri_index > -1)
//             return Hexaloop37[tetra_hex_tri_index];
//         else if (size == 3) {
//             if (tetra_hex_tri_index > -1)
//                 return Triloop37[tetra_hex_tri_index];
//             return (energy + (type>2 ? TerminalAU37 : 0));
//         }
//     // }
// #endif

    energy += mismatchH37[type][nuci1][nucj_1];

    return energy;
}

inline int v_score_hairpin_min(int i, int j, int tetra_hex_tri_index = -1) {
    int nuci = 3, nuci1 = 3, nucj_1 = 1, nucj = 2;
    int size = j-i-1;
    int type = NUC_TO_PAIR(nuci, nucj);

    int energy;

    if(size <= 30)
        energy = hairpin37[size];
    else
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));

    if(size < 3) return energy; /* should only be the case when folding alignments */
// #ifdef SPECIAL_HP
//     // if(special_hp){
//         if (size == 4 && tetra_hex_tri_index > -1)
//             return Tetraloop37[tetra_hex_tri_index];
//         else if (size == 6 && tetra_hex_tri_index > -1)
//             return Hexaloop37[tetra_hex_tri_index];
//         else if (size == 3) {
//             if (tetra_hex_tri_index > -1)
//                 return Triloop37[tetra_hex_tri_index];
//             return (energy + (type>2 ? TerminalAU37 : 0));
//         }
//     // }
// #endif

    energy += mismatchH37[type][nuci1][nucj_1];

    return energy;
}

inline int v_score_single_max(int i, int j, int p, int q){
    int type = NUC_TO_PAIR(1, 1);
    int type_2 = NUC_TO_PAIR(1, 1);
    int nucj_1 = 1, nuci1 = 1, nucp_1 = 1, nucp = 1, nucq = 1, nucq1 = 1;
    
    int n1 = p-i-1;
    int n2 = j-q-1;
    int nl, ns, u, energy;
    energy = 0;

    if (n1>n2) { nl=n1; ns=n2;}
    else {nl=n2; ns=n1;}

    if (nl == 0) {
      int nuci = 3, nucj = 4, nucp = 4, nucq = 3;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);
      return stack37[type][type_2];  /* stack */
    }

    if (ns==0) {                      /* bulge */
        energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));

    if (nl==1) {
      int nuci = 3, nucj = 4, nucp = 4, nucq = 3;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);
      energy += stack37[type][type_2];
    }
    else {
      energy += TerminalAU37;
      energy += TerminalAU37;
    }

    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1) {                   /* 1x1 loop */
        int nuci=2, nuci1=1, nucj_1=1, nucj=3, nucp=3, nucq=4;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        return int11_37[type][type_2][nuci1][nucj_1];
      }
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1) {
          int nuci=1, nucj=4, nucp=1, nucq=4, nuci1=1, nucj_1=1, nucq1=1;
          int type = NUC_TO_PAIR(nuci, nucj);
          int type_2 = NUC_TO_PAIR(nucq, nucp);
          energy = int21_37[type][type_2][nuci1][nucq1][nucj_1];
        } else {
          int nuci=1, nucj=4, nucp=1, nucq=4, nuci1=1, nucp_1=1, nucq1=1;
          int type = NUC_TO_PAIR(nuci, nucj);
          int type_2 = NUC_TO_PAIR(nucq, nucp);
          energy = int21_37[type_2][type][nucq1][nuci1][nucp_1];
        }
        return energy;
      }
      else {  /* 1xn loop */
        int nuci=1, nucj=4, nucp=1, nucq=4, nuci1=1, nucj_1=1, nucp_1=1, nucq1=1;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        int nuci=3, nucj=4, nucp=4, nucq=3, nuci1=1, nucj_1=1, nucp_1=4, nucq1=4;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1];}
      else if (nl==3){              /* 2x3 loop */
        int nuci=1, nucj=4, nucp=1, nucq=4, nuci1=1, nucj_1=1, nucp_1=1, nucq1=1;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      int nuci=1, nucj=4, nucp=1, nucq=4, nuci1=1, nucj_1=1, nucp_1=1, nucq1=1;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1];
    }
  }
  return energy;
}

inline int v_score_single_min(int i, int j, int p, int q){
    int type = NUC_TO_PAIR(1, 1);
    int type_2 = NUC_TO_PAIR(1, 1);
    int nucj_1 = 1, nuci1 = 1, nucp_1 = 1, nucp = 1, nucq = 1, nucq1 = 1;
    
    int n1 = p-i-1;
    int n2 = j-q-1;
    int nl, ns, u, energy;
    energy = 0;

    if (n1>n2) { nl=n1; ns=n2;}
    else {nl=n2; ns=n1;}

    if (nl == 0) {
      int nuci = 3, nucj = 2, nucp = 2, nucq = 3;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);
      return stack37[type][type_2];  /* stack */
    }

    if (ns==0) {                      /* bulge */
        energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));

    if (nl==1) {
      int nuci = 3, nucj = 2, nucp = 2, nucq = 3;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);
      energy += stack37[type][type_2];
    }
    else {
      // energy += TerminalAU37;
      // energy += TerminalAU37;
    }

    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1) {                   /* 1x1 loop */
        int nuci=3, nuci1=3, nucj_1=3, nucj=2, nucp=2, nucq=3;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        return int11_37[type][type_2][nuci1][nucj_1];
      }
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1) {
          int nuci=2, nucj=3, nucp=2, nucq=3, nuci1=1, nucj_1=1, nucq1=3;
          int type = NUC_TO_PAIR(nuci, nucj);
          int type_2 = NUC_TO_PAIR(nucq, nucp);
          energy = int21_37[type][type_2][nuci1][nucq1][nucj_1];
        } else {
          int nuci=3, nucj=2, nucp=3, nucq=2, nuci1=1, nucp_1=3, nucq1=3;
          int type = NUC_TO_PAIR(nuci, nucj);
          int type_2 = NUC_TO_PAIR(nucq, nucp);
          energy = int21_37[type_2][type][nucq1][nuci1][nucp_1];
        }
        return energy;
      }
      else {  /* 1xn loop */
        int nuci=2, nucj=3, nucp=2, nucq=3, nuci1=1, nucj_1=1, nucp_1=1, nucq1=1;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        int nuci=3, nucj=2, nucp=2, nucq=3, nuci1=3, nucj_1=1, nucp_1=1, nucq1=3;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1];}
      else if (nl==3){              /* 2x3 loop */
        int nuci=3, nucj=2, nucp=2, nucq=3, nuci1=3, nucj_1=1, nucp_1=1, nucq1=3;
        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      int nuci=2, nucj=3, nucp=2, nucq=3, nuci1=3, nucj_1=1, nucp_1=1, nucq1=3;
      int type = NUC_TO_PAIR(nuci, nucj);
      int type_2 = NUC_TO_PAIR(nucq, nucp);

      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1];
    }
  }
  return energy;
}

inline int v_score_M1_max(int i, int j, int k, int len, int dangle_model) {
  int nuci, nuck, nuci_1, nuck1;
  if (!dangle_model)
    nuci = 1, nuck = 4, nuci_1 = 1, nuck1 = 1;
  else
    nuci = 4, nuck = 3, nuci_1 = 1, nuck1 = 1;

  int p = i;
  int q = k;
  int tt = NUC_TO_PAIR(nuci, nuck);

  return E_MLstem(tt, nuci_1, nuck1, dangle_model);
}

inline int v_score_M1_min(int i, int j, int k, int len, int dangle_model) {
    int nuci, nuck, nuci_1, nuck1;
    if (!dangle_model)
      nuci = 2, nuck = 3, nuci_1 = 1, nuck1 =1;
    else
      nuci = 2, nuck = 3, nuci_1 = 1, nuck1 =3;

    int p = i;
    int q = k;
    int tt = NUC_TO_PAIR(nuci, nuck);

    return E_MLstem(tt, nuci_1, nuck1, dangle_model);
}

inline int v_score_multi_max(int i, int j, int len, int dangle_model) {
  int nuci, nucj, nuci1, nucj_1;
  if (!dangle_model)
    nuci = 1, nucj = 4, nuci1 = 1, nucj_1 = 1;
  else
    nuci = 3, nucj = 4, nuci1 = 1, nucj_1 = 1;

  int tt = NUC_TO_PAIR(nucj, nuci); // lhuang: closing pair in multi: reversed

    return E_MLstem(tt, nucj_1, nuci1, dangle_model) + ML_closing37;
}

inline int v_score_multi_min(int i, int j, int len, int dangle_model) {
  int nuci, nucj, nuci1, nucj_1;
  if (!dangle_model)
    nuci = 2, nucj = 3, nuci1 = 1, nucj_1 =1;
  else
    nuci = 2, nucj = 3, nuci1 = 3, nucj_1 =3;

  int tt = NUC_TO_PAIR(nucj, nuci); // lhuang: closing pair in multi: reversed

  return E_MLstem(tt, nucj_1, nuci1, dangle_model) + ML_closing37;
}

// exterior_loop
inline int v_score_external_paired_max(int i, int j, int len, int dangle_model) {
  int nuci, nucj, nuci_1, nucj1;
  if (!dangle_model)
    nuci = 1, nucj = 4, nuci_1 = 1, nucj1 =1;
  else
    nuci = 4, nucj = 3, nuci_1 = 1, nucj1 =1;

  int type = NUC_TO_PAIR(nuci, nucj);
  int energy = 0;

  if (dangle_model != 0){

      if(nuci_1 >= 0 && nucj1 >= 0){
          energy += mismatchExt37[type][nuci_1][nucj1];
      }
      else if (nuci_1 >= 0){
          energy += dangle5_37[type][nuci_1];
      }
      else if (nucj1 >= 0){
          energy += dangle3_37[type][nucj1];
      }

  }

  if(type > 2)
      energy += TerminalAU37;
  return energy;
}

inline int v_score_external_paired_min(int i, int j, int len, int dangle_model) {
  int nuci, nucj, nuci_1, nucj1;
  if (!dangle_model)
    nuci = 2, nucj = 3, nuci_1 = 1, nucj1 =1;
  else
    nuci = 2, nucj = 3, nuci_1 = 1, nucj1 =3;

  int type = NUC_TO_PAIR(nuci, nucj);
  int energy = 0;

  if (dangle_model != 0){

      if(nuci_1 >= 0 && nucj1 >= 0){
          energy += mismatchExt37[type][nuci_1][nucj1];
      }
      else if (nuci_1 >= 0){
          energy += dangle5_37[type][nuci_1];
      }
      else if (nucj1 >= 0){
          energy += dangle3_37[type][nucj1];
      }

  }

  if(type > 2)
      energy += TerminalAU37;
  return energy;
}
#endif //FASTCKY_UTILITY_V_MAX_H
