uniform vec3 light;
varying vec3 position;
varying vec3 normal;

float diffuse( vec3 N, vec3 L )
{
    return max( 0., dot( N, L ));
}

void main()
{
    float d = 1. - gl_Color.r;
    float r = (1. - d*d) * .8;
    float g = (1. - (2. * (d - .5)) * (2. * (d - .5))) * .7;
    float b = (1. - (1. - d) * (1. - d));
    vec3 color = vec3(r, g, b);

    float h = gl_Color.r;
    h = h * 10.;
    h = h - floor( h );
    h = (1. / (1. + exp(-100.*(h - .55)))) + (1. / (1. + exp(-100.*(-h + .45))));
    h = 1. - h;
    color.xyz = vec3(h, h, h) + (1. - h) * color.xyz;

    vec3 N = normalize( normal );
    vec3 L = normalize( light - position );

    gl_FragColor.rgb = diffuse(N,L) * color;
    gl_FragColor.a = 1.0;
}
