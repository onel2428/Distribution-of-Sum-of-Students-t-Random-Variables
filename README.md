# Distribution-of-Sum-of-Students-t-Random-Variables
Two accurate methods for fitting the distribution linear combinations of Student's t random variables with more than 2 degrees of freedom to that of a Scaled Student's t random variable. 

The first method is based on the matching of Seconds and Absolute Moments and can be implemented by using the function 
[v, s] = TFitting_KRVs_AbsoluteMoment(v,s)

The second method is based on the matching of Seconds moments and Characteristic Functions and can be implemented by using the function 
[vf, sf, rf] = TFitting_KRVs_CharacteristicFunction(v,s,r)

See functions' description for help on how to use them. The following functions are auxiliar and used for those above
M = SecondMoment(nu,sigma)       --> second moment of the sum of t-RVs with vectors nu and sigma specifying their number of degrees of freedom and scale 
M = AbsoluteMoment_2RVs(nu,sigma)--> second moment of the sum of two t-RVs with two-element vectors nu and sigma specifying their number of degrees of freedom and scale  

There is also a function that can be called to assess the fitting accuracy in terms of the Bhattacharyya distance:
divBha = Bhattacharyya_distance(t,s,v)
