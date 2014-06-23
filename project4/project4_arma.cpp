#include <cstdio>
#include "../utils.hpp"

using namespace std;

#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

void mp2_noddy(arma::vec& TEI_MO, const arma::vec& TEI_AO, const arma::mat& C)
{
  int NBasis = C.n_rows;

  int i, j, k, l, ijkl;
  int p, q, r, s, pq, rs, pqrs;

  for (i = 0, ijkl = 0; i < NBasis; i++) {
    for (j = 0; j <= i; j++) {
      for (k = 0; k <= i; k++) {
	for (l = 0; l <= (i==k ? j : k); l++, ijkl++) {
	  for (p = 0; p < NBasis; p++) {
	    for (q = 0; q < NBasis; q++) {
	      pq = INDEX(p,q);
	      for (r = 0; r < NBasis; r++) {
		for (s = 0; s < NBasis; s++) {
		  rs = INDEX(r,s);
		  pqrs = INDEX(pq,rs);
		  TEI_MO(ijkl) += C(p,i) * C(q,j) * C(r,k) * C(s,l) * TEI_AO(pqrs);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return;
}

void mp2_smart(arma::vec& TEI_MO, const arma::vec& TEI_AO, const arma::mat& C)
{
  int NBasis = C.n_rows;
  int M = NBasis*(NBasis+1)/2;
  
  arma::mat tmp = arma::mat(M, M, arma::fill::zeros);
  arma::mat X = arma::mat(NBasis, NBasis, arma::fill::zeros);

  int i, j, k, l, ij, kl, ijkl, klij;

  // for(i=0,ij=0; i < nao; i++)
  //   for(j=0; j <= i; j++,ij++) {
  //     for(k=0,kl=0; k < nao; k++)
  // 	for(l=0; l <= k; l++,kl++) {
  // 	  ijkl = INDEX(ij,kl);
  // 	  X[k][l] = X[l][k] = TEI[ijkl];
  // 	}
  //     zero_matrix(Y, nao, nao);
  //     mmult(C, 1, X, 0, Y, nao, nao, nao);
  //     zero_matrix(X, nao, nao);
  //     mmult(Y, 0, C, 0, X, nao, nao, nao);
  //     for(k=0, kl=0; k < nao; k++)
  // 	for(l=0; l <= k; l++, kl++)
  // 	  TMP[kl][ij] = X[k][l];
  //   }

  // zero_array(TEI, (nao*(nao+1)/2)*((nao*(nao+1)/2)+1)/2);

  // for(k=0,kl=0; k < nao; k++)
  //   for(l=0; l <= k; l++,kl++) {
  //     zero_matrix(X, nao, nao);
  //     zero_matrix(Y, nao, nao);
  //     for(i=0,ij=0; i < nao; i++)
  // 	for(j=0; j <=i; j++,ij++)
  // 	  X[i][j] = X[j][i] = TMP[kl][ij];
  //     zero_matrix(Y, nao, nao);
  //     mmult(C, 1, X, 0, Y, nao, nao, nao);
  //     zero_matrix(X, nao, nao);
  //     mmult(Y, 0, C, 0, X, nao, nao, nao);
  //     for(i=0, ij=0; i < nao; i++)
  // 	for(j=0; j <= i; j++,ij++) {
  // 	  klij = INDEX(kl,ij);
  // 	  TEI[klij] = X[i][j];
  // 	}
  //   }

  for (i = 0, ij = 0; i < NBasis; i++)
    for (j = 0; j <= i; j++, ij++) {
      for (k = 0, kl = 0; k < NBasis; k++)
	for (l = 0; l <= k; l++, kl++)
	  X(k,l) = X(l,k) = TEI_AO(INDEX(ij,kl));
      X = C.t() * X * C;
      for (k = 0, kl = 0; k < NBasis; k++)
	for (l = 0; l <= k; l++, kl++)
	  tmp(kl,ij) = X(k,l);
    }
  
  TEI_MO.zeros();

  for (k = 0, kl = 0; k < NBasis; k++)
    for (l = 0; l <= k; l++, kl++) {
      for (i = 0, ij = 0; i < NBasis; i++)
	for (j = 0; j <= i; j++, ij++)
	  X(i,j) = X(j,i) = tmp(kl,ij);
      X = C.t() * X * C;
      for (i = 0, ij = 0; i < NBasis; i++)
	for (j = 0; j <= i; j++, ij++)
	  TEI_MO(INDEX(kl,ij)) = X(i,j);
    }

  return;
}

double calc_MP2_energy(const arma::vec& TEI_MO, const arma::vec& E)
{
  double E_MP2 = 0.0;

  int NBasis = E.n_rows;
  int NOcc = 5;
  int NVirt = NBasis - NOcc;

  // for (int i = 0; i < NOcc; i++)
  //   for (int j = 0; j < NOcc; j++)
  //     for (int a = NOcc; a < NVirt; a++)
  // 	for (int b = NOcc; b < NVirt; b++)
  // 	  E_MP2 += (TEI_MO(idx4(i,a,j,b)) * (2*TEI_MO(idx4(i,a,j,b)) - TEI_MO(idx4(i,b,j,a))))/(E(i) + E(j) - E(a) - E(b));

  int i, j, a, b, ia, ja, ib, jb, iajb, ibja;

  for (i = 0; i < NOcc; i++) {
    for (a = NOcc; a < NBasis; a++) {
      ia = INDEX(i,a);
      for (j = 0; j < NOcc; j++) {
	ja = INDEX(j,a);
	for (b = NOcc; b < NBasis; b++) {
	  ib = INDEX(i,b);
	  jb = INDEX(j,b);
	  iajb = INDEX(ia,jb);
	  ibja = INDEX(ib,ja);
	  E_MP2 += TEI_MO(iajb) * (2*TEI_MO(iajb) - TEI_MO(ibja))/(E(i) + E(j) - E(a) - E(b));
	}
      }
    }
  }

  return E_MP2;
}

int main()
{
  // read in the TEIs and MO coefficients/energies previously saved to disk.
  arma::vec TEI_AO;
  bool status_TEI_AO = TEI_AO.load("ERI.mat");
  arma::mat C;
  bool status_C = C.load("C.mat");
  arma::mat F_MO;
  bool status_F_MO = F_MO.load("F_MO.mat");

  if (status_TEI_AO == true) printf("TEI loading worked.\n");
  else printf("TEI loading failed.\n");
  if (status_C == true) printf("C loading worked.\n");
  else printf("C loading failed.\n");
  if (status_F_MO == true) printf("F_MO loading worked.\n");
  else printf("F_MO loading failed.\n");

  arma::vec TEI_MO = arma::vec(TEI_AO.n_elem, arma::fill::zeros);
  // mp2_noddy(TEI_MO, TEI_AO, C);
  mp2_smart(TEI_MO, TEI_AO, C);

  // TEI_AO.print("TEI_AO:");
  // TEI_MO.print("TEI_MO:");

  arma::vec E = F_MO.diag();
  // printf("HF MO energies:\n"); print_arma_vec(E, 3);
  // this is tr(PF), not tr(P(H+F))!
  // double E_HF = arma::sum(E);
  // printf(" E_HF: %20.12f\n", E_HF);

  double E_MP2 = calc_MP2_energy(TEI_MO, E);
  printf("E_MP2: %20.12f\n", E_MP2);

  return 0;
}
