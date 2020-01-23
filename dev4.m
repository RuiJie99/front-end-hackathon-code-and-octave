
function [xi, yi, zi, face]=Devoir4(nout,nin,poso)

  %nous ferons pour chaque rayon le m�me processus. On regarde s'il frappe
  %l'objet, et si oui, quelle image virtuelle cela contribuera a former. 
  %Nous allons tester une quantit� n suffisant de rayon allant des angles
  %polaires a azitumales. Pour l'instant nous avons mis i=10000. Et nous
  %allons aller de -pi/deux � pi/deux. 
  %Nous allons mettre les positions d'image virtuelles x,y,z dans les
  %vecteurs suivants:
  xi=[];
  yi=[];
  zi=[];
  %La face touch�e par les rayons sera inscrite dans le vecteur suivant:
  face=[];
  dist_parcourue = 0;

  R_ellipse=3;

  %calcul distance entre surface et foyers
  r_foyer1=[4;4;2.514];
  r_foyer2=[4;4;19.485];
  r_top = [4; 4; 20];
  dist_foyers = norm(r_foyer1 - r_top) + norm(r_foyer2 - r_top);

  p_bloc = [4; 4; 11];
  % Angles vertical de l'ellipse 
  [ angle_fin_phi,angle_debut_phi] = CalculerAnglesCercle(poso);
  [angle_fin_theta,angle_debut_theta] = CalculerAnglesEllipse(poso);
  n=100;
  dist = 0.1;

  variationAnglePhi=(angle_fin_phi-angle_debut_phi)/n;
  variationAngleTheta=(angle_fin_theta-angle_debut_theta)/n;
  theta = angle_debut_theta;

      while( theta <= angle_fin_theta)  
          phi = angle_debut_phi;
          while (phi <= angle_fin_phi)

                dir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
                RayonAManqueEllipse=0;
                RayonAToucheEllipse=0;
                r_rayon=poso;

                while(RayonAToucheEllipse==0 && RayonAManqueEllipse==0)
                %avancer le rayon
                  r_rayon_previous = r_rayon;
                  r_rayon = AvancerRayon(r_rayon, dir, dist);
                  % Le rayon a t-il touch� l'ellipse ?
                  RayonAToucheEllipse = (norm(r_rayon - r_foyer1) + norm(r_rayon - r_foyer2)) <= dist_foyers;
                  % Depasse le haut et continue de monter ou depasse le bas et continue de descendre
                  RayonAManqueEllipseZ = (r_rayon(3)<2 && r_rayon(3) <= r_rayon_previous(3)) || (r_rayon(3) > 20 && r_rayon(3) >= r_rayon_previous(3)) && RayonAToucheEllipse == 0;
                  % S'eloigne de la surface du cylindre
                  RayonAManqueEllipseXY = (norm(r_rayon - r_foyer1) + norm(r_rayon - r_foyer2)) > (norm(r_rayon_previous - r_foyer1) + norm(r_rayon_previous - r_foyer1)) && RayonAToucheEllipse == 0;
                  RayonAManqueEllipse = RayonAManqueEllipseXY || RayonAManqueEllipseZ;
                end
                dir_refrac = [];
                %Si le rayon n'a pas �t� r�fl�chi, on continue et on regarde si le
                %rayon fini par toucher le bloc de m�tal
                rayonAToucheLeBlocMetal=0;
                compteurRefractionInterne=0;
                rayonEstSortiEllipse=0;
                PositionAvantEntrerDansEllipse=r_rayon;
                if(RayonAToucheEllipse==1)
                    %Angle entrant

                    angleEntrant = acos(dot(dir, N)/(norm(dir)*norm(N)));
                    angleRefracte = asin(nin/nout * sin(angleEntrant));
                    si2 = nout/nin * sin(angleEntrant);
                    %mettre le nouvel angle post-ellipse r�fract�
                    if (abs(si2) < 1)
                      while(rayonAToucheLeBlocMetal==0 && rayonEstSortiEllipse==0)
                        rayonAToucheLeBlocMetal=DetectionToucheBlocMetal(r_rayon); %rapporte 0 si n'a pas touch� ou sinon le num�ro de la face  touch�e
                        d_rayon=norm(r_rayon-r_foyer1)+norm(r_rayon-r_foyer2);
                        d_r = norm([r_rayon(1) - p_bloc(1); r_rayon(2)-p_bloc(2)]);
                        r = (4 + sqrt(9*(1 - (((r_rayon(3) - 11)^2)/(9^2))))) - p_bloc(1);
                        if(d_r >= r)%Cas sp�cial si le n de l'ellipse est plus grand qu'� l'ext�rieur, qui permet les r�flections internes
                            %On regarde si le rayon est arriv� � un point ou il y a
                            %r�flexion interne
                            %updater l'angle selon la nouvelle reflexion
                          j = cross(dir_refrac, N)/norm(cross(dir_refrac, N));
                          k = cross(N, j);
                          angleEntrant = acos(dot(dir_refrac, N)/(norm(dir_refrac)*norm(N)));
                          dir_refrac = -N*sqrt(1-(nin/nout*dot(dir_refrac, k))^2) + k*(nin/nout*dot(dir_refrac, k));
                          si2 = nin/nout * sin(angleEntrant);
                          compteurRefractionInterne=compteurRefractionInterne+1;
                        end
                        rayonEstSortiEllipse= abs(si2) < 1 || compteurRefractionInterne==100; %%on rajoute 1 de distance pour �tre 100% sur qu'il est sorti et qu'il a aurait eu le temps de faire la r�flexion interne si c'est le cas
                        r_rayon = AvancerRayon(r_rayon, dir_refrac, dist);
                        dist_parcourue =dist_parcourue+ dist * dir_refrac/norm(dir_refrac);
                        si2 = 1;
                      end
                    end
                end
                if (rayonAToucheLeBlocMetal >=1)
                    face=[face; rayonAToucheLeBlocMetal];
                    posImageVirtuelle=AvancerRayon(PositionAvantEntrerDansEllipse,dir,norm(dist_parcourue));
                    xi=[xi; posImageVirtuelle(1)];
                    yi=[yi; posImageVirtuelle(2)];
                    zi=[zi; posImageVirtuelle(3)];
                end

              phi = phi + variationAnglePhi;
          end
          theta = theta + variationAngleTheta; 
      end    
end

% ------------------------------------ %

% RJ
function new_dir = ReflechirRayon(pos, dir)
  r_rayon = pos;
  % Normal au plan
  N = [(r_rayon(1)-4)/(9);(r_rayon(2)-4)/(9);(r_rayon(3)-11)/(81)]/norm([(r_rayon(1)-4)/(9);(r_rayon(2)-4)/(9);(r_rayon(3)-11)/(81)]);
  disp(N);
  % Normal � la direction et dans le plan
  j = cross(dir, N)/norm(cross(dir, N));
  disp(j);
  % vecteur k dans le plan et perpendiculaire au vecteur j
  k = cross(N, j);
  disp(k);
  % sinus direction de transmission
  sr = dot(dir/norm(dir),k);
  disp(sr);
  % new direction
  new_dir = N*sqrt(1 - sr^2) + k*sr;
end

% RJ
function new_dir = RefracterRayon(pos, dir, ni, nt)
  r_rayon = pos;
  if (ni == nt)
      new_dir = dir;
  else
      % Normal au plan
      N = [(r_rayon(1)-4)/(9);(r_rayon(2)-4)/(9);(r_rayon(3)-11)/(81)]/norm([(r_rayon(1)-4)/(9);(r_rayon(2)-4)/(9);(r_rayon(3)-11)/(81)]);
      % Normal � la direction et dans le plan
      j = cross(dir, N)/norm(cross(dir, N));
      % vecteur k dans le plan et perpendiculaire au vecteur j
      k = cross(N, j);
      % sinus direction de transmission
      st = ni/nt * dot(dir,k);
      % new direction
      new_dir = - N*sqrt(1 - st^2) + k*st;
  end

end


% RJ
function [angle_max, angle_min]=CalculerAnglesEllipse(posObs) 
  % Projeter la position de l'observateur et celle de l'ellisoide sur le plan XZ
  % Tous les points sont deplacer de 4 vers la gauche et de 11 vers le bas
  % Ainsi l'ellipsoide est centrer a (0, 0)
  rx_rz = [posObs(1)-4; posObs(3)-11];

  a = -(9*rx_rz(2)^2+81*rx_rz(1)^2);
  b = 1458*rx_rz(1);
  c = (81*rx_rz(2)^2-6561);
  %Positions des x
  x_1 = (-b + sqrt(power(b,2)-4*a*c))/(2*a);
  x_2 = (-b - sqrt(power(b,2)-4*a*c))/(2*a);

  % Les angles de la moitie en haut : y = sqrt(81-9*x_1^2)
  angle_1_haut = atan(((rx_rz(2)-sqrt(81-9*x_1^2))/((rx_rz(1)-x_1))));
  angle_2_haut = atan(((rx_rz(2)-sqrt(81-9*x_2^2))/((rx_rz(1)-x_2))));

  % Les angles de la moitie en bas : y = -sqrt(81-9*x_1^2)
  angle_1_bas =  atan((rx_rz(2)+sqrt(81-9*x_1^2))/((rx_rz(1)-x_1)));
  angle_2_bas =  atan((rx_rz(2)+sqrt(81-9*x_2^2))/((rx_rz(1)-x_2)));

  angle_max = max([angle_1_haut, angle_2_haut, angle_1_bas, angle_2_bas]);
  angle_min = min([angle_1_haut, angle_2_haut, angle_1_bas, angle_2_bas]);
end

% RJ
function [angle_max, angle_min]=CalculerAnglesCercle(posObs)
  %Position en x et en y
  rx_ry = [posObs(1)-4; posObs(2)-4];
  % Les valeurs de ax^2 + bx + c = 0 comme calcul� dans le word
  a = -(rx_ry(1)^2 + rx_ry(2)^2);
  b = 18*rx_ry(1);
  c = (9*rx_ry(2)^2-81);
  %Positions des x
  x_1 = (-b + sqrt(power(b,2)-4*a*c))/(2*a);
  x_2 = (-b - sqrt(power(b,2)-4*a*c))/(2*a);

  % Les angles de la moitie en haut : y = sqrt(81-9*x_1^2)
  angle_1_haut = atan((rx_ry(2)-sqrt(9-x_1^2))/((rx_ry(1)-x_1)));
  angle_2_haut = atan((rx_ry(2)-sqrt(9-x_2^2))/((rx_ry(1)-x_2)));

  % Les angles de la moitie en bas : y = -sqrt(81-9*x_1^2)
  angle_1_bas = atan((rx_ry(2)+sqrt(9-x_1^2))/((rx_ry(1)-x_1)));
  angle_2_bas = atan((rx_ry(2)+sqrt(9-x_2^2))/((rx_ry(1)-x_2)));

  angle_max = max([angle_1_haut, angle_2_haut, angle_1_bas, angle_2_bas]);
  angle_min = min([angle_1_haut, angle_2_haut, angle_1_bas, angle_2_bas]);
end

% ---------------------------------- %

function resultat=DetectionToucheBlocMetal(r_rayon)
  resultat = 0;
  pos_bloc = [3.5; 4; 14.5];
  h_x_b = 0.5;
  h_y_b = 1;
  h_z_b = 2.5;
  f1 = [3; 0; 0]; %rouge
  f2 = [4; 0; 0]; %cian
  f3 = [0; 3; 0]; %vert
  f4 = [0; 5; 0]; %jaune
  f5 = [0; 0; 12]; %bleu
  f6 = [0; 0; 17]; %magenta
  n1 = f1/norm(f1);
  n2 = f2/norm(f2);
  n3 = f3/norm(f3);
  n4 = f4/norm(f4);
  n5 = f5/norm(f5);
  n6 = f6/norm(f6);
  dn1 = (dot(n1, (f1 - transpose(r_rayon)))/norm(n1));
  dn2 = (dot(n2, (f2 - transpose(r_rayon)))/norm(n2));
  dn3 = (dot(n3, (f3 - transpose(r_rayon)))/norm(n3));
  dn4 = (dot(n4, (f4 - transpose(r_rayon)))/norm(n4));
  dn5 = (dot(n5, (f5 - transpose(r_rayon)))/norm(n5));
  dn6 = (dot(n6, (f6 - transpose(r_rayon)))/norm(n6));
  dn = [dn1; dn2; dn3; dn4; dn5; dn6];
  dist_f = min([abs(dn1) abs(dn2) abs(dn3) abs(dn4) abs(dn5) abs(dn6)]);
  face = 0;
  for (i = 1:6)
    if (abs(dn(i)) == dist_f)  
        face = i; %%Identifie la face
        break;
    end
  end
  if (dist_f < 0.1)
      if ((((face == 1 && dn1 >= 0) || (face == 2 && dn2 <= 0)) &&(r_rayon(2) < (pos_bloc(2) + h_y_b)) && (r_rayon(2) > (pos_bloc(2) - h_y_b)) &&(r_rayon(3) < (pos_bloc(3) + h_z_b)) && (r_rayon(3) > (pos_bloc(3) - h_z_b))) ||(((face == 3 && dn3 >= 0) || (face == 4 && dn4 <= 0))&& (r_rayon(1) < (pos_bloc(1) + h_x_b)) && (r_rayon(1) > (pos_bloc(1) - h_x_b))&&(r_rayon(3) < (pos_bloc(3) + h_z_b)) && (r_rayon(3) > (pos_bloc(3) - h_z_b)))||(((face == 5 && dn5 >= 0) || (face == 6 && dn6 <= 0))&&(r_rayon(2) < (pos_bloc(2) + h_y_b)) && (r_rayon(2) > (pos_bloc(2) - h_y_b))&&(r_rayon(1) < (pos_bloc(1) + h_x_b)) && (r_rayon(1) > (pos_bloc(1) - h_x_b))))
          resultat = face;
          fprintf("face %d\n", face);
      end
  end
end


function new_pos = AvancerRayon(pos, direction, distance)
  unit = direction/norm(direction);
  new_pos = pos + transpose(unit)*distance;
end