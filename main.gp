/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
gen = ffgen(q,'a);
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);




/* ****************************************************************************** */
/* ******************************* Explications : ******************************* */
/* ****************************************************************************** */


/*
 * L'interlocuteur d'Alice, nommons-le Buzz, a chiffré son message,
   sous la forme d'un code BCH.
 * Il a ensuite introduit des erreurs afin de brouiller les pistes.

 * On connaît la distance prévue delta = 2t+1.
 
 * Il nous faut décoder le message BCH en corrigeant les erreurs.
 * Pour ce faire, on va utiliser
   * le polynôme syndrome S(x) = Sum_{l=0}^{2t-1} y(alpha^{a+l})*x^l,
   * le polynôme localisateur d'erreur E(x) = Prod_{i \in T} (1- alpha^i *x)
   * et le polynôme R(x) = sum_{i \in T} e_i * alpha^{ia} Prod_{j \in (T\i)} (1-alpha^j * x).
 * Les deux derniers polynômes sont obtenus comme solutions de
   l'équation de Padé S(x)E(x)=R(x) [x^{2t}].

 * Les racines du polynôme E permettent de retrouver les indices des erreurs.
 
*/




/* ****************************************************************************** */
/* ******************************** Fonctions : ********************************* */
/* ****************************************************************************** */

/* Pour obtenir y(alpha^{a+l}), on substitue x' dans y par alpha^(a+l). */

get_S(y, alpha, a, t)={ return (sum(l=0, 2*t-1, subst(y, 'x, alpha^(a+l)) * 'x^l) ); };



/* Comme décrit ci-dessus, on cherche à obtenir la solution de l'équation de Padé
   S(x)E(x)=R(x) [x^{2t}].
   Pour ce faire, on utilise la fonction bestapprPade(x, {B}).
 * bestapprPade(x, {B}): returns a rational function approximation to x. This 
   function applies to series, polmods, and rational functions of course. 
   Otherwise it applies recursively to all components.
*/
get_Pade_approx(y, alpha, a, t)={
        my(synd, pade);
        synd = get_S(y, alpha, a, t);
        pade = bestapprPade(Mod(synd, x^(2*t)));
        return ([numerator(pade), denominator(pade)]);
};


/*
 * Une fois E récupéré, ses racines (2t-1 puissances de alpha consécutives)
   permettent de retrouver les indices des erreurs :
 * E(a^{-k})=0 indique qu'il y a une erreur en indice i : on la trouve en calculant R(a^{-k})/E(a^{-k}).

 * Cependant, on ne sait pas quelle est la première puissance de alpha telle qu'elle soit une racine.
 * Pour les trouver, on parcourt l'ensemble des racines primitives alpha jusqu'à ce qu'une d'entre elles
   vérfie les conditions (y(alpha^(a+l)) = 0 pour tout mot y du code et tout l dans {1,...,2t}).

 * On note dans une liste les erreurs détectées pour la valeur de alpha testée.
 * Enfin, on suppose qu'un nombre suffisamment grand d'erreurs représente un message :
   si la taille minimale est respectée, on renvoie le message.
*/

get_e_k(y, min_taille_message)={
        my(R,E, e_k, alpha, val, R_E_);
	for(i = 0, q,
	     alpha = ffprimroot(gen);
	     for(a = 0, q-1,
	       [R,E] = get_Pade_approx(y, alpha, a, t);
	       e_k = List();
	       for(k = 0, q-2, R_E_ = subst(E,'x,alpha^(-k));
	       if(R_E_ == 0, val = subst((R/deriv(E))*(x^(a-1)), 'x , alpha^(-k));
	         listput(e_k, fqx2int(val, gen)));
	       );
	       if(#e_k > min_taille_message,  return (Strchr(Vec(e_k))));
	        );
	 );
};




/* ****************************************************************************** */
/* ******************************* Application : ******************************** */
/* ****************************************************************************** */


y = int2fqx(m, gen);
e_ k = get_e_k(y, 6);
print(e_k);