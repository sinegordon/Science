
# define RNGInit unsigned int x1, y1, z1, w1; x1 = seed[idx]; y1 = 362436069; z1 = 521288629; w1 = 88675123
# define rand_uniform_float random_uniform_float(&x1, &y1, &z1, &w1)
# define rand_uniform_int random_uniform_int(&x1, &y1, &z1, &w1)

inline float random_uniform_float(unsigned int * x1, unsigned int * y1, unsigned int * z1, unsigned int * w1)
{
    unsigned int t1 = ((*x1)^((*x1)<<11));
    (*x1) = (*y1);
    (*y1) = (*z1);
    (*z1) = (*w1); 
    (*w1) = ((*w1)^((*w1)>>19))^(t1^(t1>>8)); 
    return ((float)(*w1))/(UINT_MAX);
}

inline unsigned int random_uniform_int(unsigned int * x1, unsigned int * y1, unsigned int * z1, unsigned int * w1)
{
    unsigned int t1 = ((*x1)^((*x1)<<11));
    (*x1) = (*y1);
    (*y1) = (*z1);
    (*z1) = (*w1); 
    (*w1) = ((*w1)^((*w1)>>19))^(t1^(t1>>8)); 
    return *w1;
}

// Вспомогательные функции для проверки пересечения двух отрезков
inline float max2(float a, float b)
{
    return a > b? a: b;
}

inline float min2(float a, float b)
{
    return a < b? a: b;
}

inline float area_sign (float ax, float ay, float bx, float by, float cx, float cy) {
	return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax) > 0 ? 1: -1;
}
 
inline bool intersect_1 (float a, float b, float c, float d) {
    float t = 0;
	if (a > b)  {t = b; a = b; b = t;};
	if (c > d)  {t = d; c = d; d = t;};
	return max2(a,c) <= min2(b,d);
}

inline bool point_in_poly1(__global float * lattice, unsigned int count_vertex, float xp, float yp)
        {
            unsigned int i, j;
            bool c = false;
            float xi, yi, xj, yj;
            for (i = 0, j = count_vertex - 1; i < count_vertex; j = i++)
            {
                xi = lattice[2*i];
                yi = lattice[2*i + 1];
                xj = lattice[2*j];
                yj = lattice[2*j + 1];
                if ((((yi <= yp) && (yp < yj)) ||
                    ((yj <= yp) && (yp < yi))) &&
                    (xp < (xj - xi) * (yp - yi) / (yj - yi) + xi))
                    c = !c;
            }
            return c;
        }
/*
inline bool point_in_poly(__global float * lattice, unsigned int count_vertex, float xp, float yp)
{
    bool result = false;
    int j = count_vertex - 1;
    float xi, yi, xj, yj;
    for (unsigned int i = 0; i < count_vertex; i++) {
        xi = lattice[2*i];
        yi = lattice[2*i + 1];
        xj = lattice[2*j];
        yj = lattice[2*j + 1];
        if ( (yi < yp && yj >= yp || yj < yp && yi >= yp ) && ( xi + (yp - yi) / (yj - yi) * (xj - xi) < xp) )
        {
            result = !result;
        };
        j = i;
    }
}
*/
// Проверка пересечения двух отрезков
inline bool eges_intersect (float ax, float ay, float bx, float by, float cx, float cy, float dx, float dy) {
    float eps = 1e-7;
	return intersect_1 (ax, bx, cx, dx)
		&& intersect_1 (ay, by, cy, dy)
		&& area_sign(ax, ay, bx, by, cx, cy) * area_sign(ax, ay, bx, by, dx, dy) <= eps
		&& area_sign(cx, cy, dx, dy, ax, ay) * area_sign(cx, cy, dx, dy, bx, by) <= eps;
}

__kernel void Step(
        __global float * lattice, 
        __global float * probe_lattice,
        __global float * centers,
        __global float * params,
        __global int * accept,
        __global int * seed
        /*
        int size, 
        int count_poly, 
        float xmin, 
        float xmax, 
        float ymin, 
        float ymax
        */
        )
{
    // Инициализируем генератор псевдослучайных чисел
    unsigned int idx = get_global_id(0);
	RNGInit;
    seed[idx] = x1;
    // Копируем параметры задачи в локальные переменные
    unsigned int size = (unsigned int) params[0];
    unsigned int count_poly = (unsigned int) params[1]; 
    float xmin = params[2];
    float xmax = params[3];
    float ymin = params[4];
    float ymax = params[5];
    // Определяем номер сдвигаемого полигона
    idx = rand_uniform_int % count_poly;
    bool flag;
    // Координат на полигон в массиве
	unsigned int coords_per_poly = size / count_poly;
	// Максимальные смещения
	float max_translation = 0.1;
	float max_angle = 0.1;
    // Вспомогательные переменные
	int local_accept = 0;
    float xc, yc, xc1, yc1;
    float x, y;
    float dist, angle;
    float xcl = 0, ycl = 0;
    float dx, dy;
    
    // Двигаем полигон с индексом k (с учетом периодичности границ)
    // Сохраняем старые координаты в probe_lattice
    for(unsigned int i = 0; i < coords_per_poly; i+=2)
    {
        probe_lattice[idx*coords_per_poly + i] = lattice[idx*coords_per_poly + i];
        probe_lattice[idx*coords_per_poly + i + 1] = lattice[idx*coords_per_poly + i + 1];
    }
    xc = centers[2*idx];
    yc = centers[2*idx + 1];
    // Вращаем
    if( rand_uniform_float > 0.5)
    {
        angle = max_angle*(1 - 2*rand_uniform_float);
        for(unsigned int i = 0; i < coords_per_poly; i+=2)
        {
            x = lattice[idx*coords_per_poly + i] - xc;
            y = lattice[idx*coords_per_poly + i + 1] - yc;
            lattice[idx*coords_per_poly + i] = x*cos(angle) - y*sin(angle) + xc;
            lattice[idx*coords_per_poly + i + 1] = x*sin(angle) + y*cos(angle) + yc;
        }
    }
    // Параллельно переносим
    else{
        dx = max_translation*(1 - 2*rand_uniform_float);
        dy = max_translation*(1 - 2*rand_uniform_float);
        if (xc + dx > xmax)
            dx = dx - (xmax - xmin);
        else if (xc + dx < xmin)
            dx = dx + (xmax - xmin);
        if (yc + dy > ymax)
            dy = dy - (ymax - ymin);
        else if (yc + dy < ymin)
            dy = dy + (ymax - ymin);
        // Параллельный перенос
        for(unsigned int i = 0; i < coords_per_poly; i+=2)
        {
            lattice[idx*coords_per_poly + i] += dx;
            lattice[idx*coords_per_poly + i + 1] += dy;
        };
        centers[2*idx] += dx;
        centers[2*idx + 1] += dy;
        xc = centers[2*idx];
        yc = centers[2*idx + 1];
    }
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    // Проверяем, что полигон с индексом idx (трансформированный) ни с кем не пересекается
    flag = false;
    dist = 0;
    for(unsigned int i = 0; i < count_poly; i++)
    {
        xc1 = centers[2*i];
        yc1 = centers[2*i + 1];
        dist = (xc - xc1)*(xc - xc1) + (yc - yc1)*(yc - yc1);
        if ( i == idx || dist > 4 )
        {
            continue;
        }
        else
        {
            flag = false;
            // Проверка по пересечению ребер
            /*
            for(unsigned int l = 0; l < coords_per_poly - 2; l+=2)
            {
                for(unsigned int k = 0; k < coords_per_poly - 2; k+=2)
                {
                    flag = eges_intersect(  lattice[idx*coords_per_poly + l], 
                                            lattice[idx*coords_per_poly + l + 1],
                                            lattice[idx*coords_per_poly + l + 2],
                                            lattice[idx*coords_per_poly + l + 3],
                                            lattice[i*coords_per_poly + k],
                                            lattice[i*coords_per_poly + k + 1],
                                            lattice[i*coords_per_poly + k + 2],
                                            lattice[i*coords_per_poly + k + 3]);
                    if ( flag == true )
                        break;
                }
                if ( flag == true )
                    break;
            }
            if ( flag == true )
                break;
            */    
            // Проверка по вхождению вершин
            for(unsigned int l = 0; l < coords_per_poly; l+=2)
            {
                flag = point_in_poly1(lattice + idx*coords_per_poly, coords_per_poly/2, lattice[i*coords_per_poly + l], lattice[i*coords_per_poly + l + 1]);
                if ( flag == true )
                    break;
                flag = point_in_poly1(lattice + i*coords_per_poly, coords_per_poly/2, lattice[idx*coords_per_poly + l], lattice[idx*coords_per_poly + l + 1]);
                 if ( flag == true )
                    break;
            }
            if ( flag == true )
                break;
        };
        if ( flag == true )
            break;
    };
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    // Если пересекается возвращаем все на место из сохраненных в probe_lattice координат
    if (flag == true)
    {
        for(unsigned int i = 0; i < coords_per_poly; i++ )
            lattice[idx*coords_per_poly + i] = probe_lattice[idx*coords_per_poly + i];
        xcl = 0;
        ycl = 0;
        for(unsigned int i = 0; i < coords_per_poly; i+=2)
        {
            xcl += lattice[idx*coords_per_poly + i];
            ycl += lattice[idx*coords_per_poly + i + 1];
        }
        xcl /= coords_per_poly;
        ycl /= coords_per_poly;
        xc = 2*xcl;
        yc = 2*ycl;
        centers[2*idx] = xc;
        centers[2*idx + 1] = yc;
    }
    else
    {
        local_accept += 1;
    };
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    accept[idx] = local_accept;
    seed[idx] = x1;
};