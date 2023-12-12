var documenterSearchIndex = {"docs":
[{"location":"transform/#Transformations","page":"Trasformations","title":"Transformations","text":"","category":"section"},{"location":"transform/#Classical-Orbital-Elements-(coe)","page":"Trasformations","title":"Classical Orbital Elements (coe)","text":"","category":"section"},{"location":"transform/","page":"Trasformations","title":"Trasformations","text":"convert6_cart_to_coe\nconvert6_coe_to_cart\n\n∂convert6_cart_to_coe\n∂convert6_coe_to_cart","category":"page"},{"location":"transform/#AstroRepresentations.convert6_cart_to_coe","page":"Trasformations","title":"AstroRepresentations.convert6_cart_to_coe","text":"convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert cartesian state representation into classical orbital elements.\n\nInputs\n\nsv – state vector – L, T\nμ – center gravitational parameter  – L³/T²\n\nOutput\n\nClassical orbital elements representation of the state as a SVector{6}. \n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications. Vol. 12.  Springer Science & Business Media, 2001.\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.convert6_coe_to_cart","page":"Trasformations","title":"AstroRepresentations.convert6_coe_to_cart","text":"convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert classical orbital elements state vector to cartesian state.\n\nInputs\n\nsv – Keplerian elements – L, rad\nμ – Center's gravitational parameter  – L³/T²\n\nOutput\n\nCartesian representation of the state as a SVector{6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications. Vol. 12.  Springer Science & Business Media, 2001.\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.∂convert6_cart_to_coe","page":"Trasformations","title":"AstroRepresentations.∂convert6_cart_to_coe","text":"∂convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert cartesian state representation into classical orbital elements. Compute also the full jacobian of the elements wrt the cartesian state.\n\nInputs\n\nsv – state vector – L, T\nμ – center gravitational parameter  – L³/T²\n\nOutput\n\nClassical Orbital Elements representation of the state as a SVector{6} and its jacobian as  a SMatrix{6, 6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications. Vol. 12.  Springer Science & Business Media, 2001.\nPasquale, A. - Multiple Shooting Optimiser (MSO). Technical Note 0001, 2022.\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.∂convert6_coe_to_cart","page":"Trasformations","title":"AstroRepresentations.∂convert6_coe_to_cart","text":"∂convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert classical orbital elements state vector to cartesian state. Compute also the full jacobian of the cartesian states wrt the elements.\n\nInputs\n\nsv – cartesian state representation – L, rad\nμ – Center's gravitational parameter  – L³/T²\n\nOutput\n\nCartesian representation of the state as a SVector{6} and its jacobian as a SMatrix{6, 6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications. Vol. 12.  Springer Science & Business Media, 2001.\nPasquale, A. - Multiple Shooting Optimiser (MSO). Technical Note 0001, 2022.\n\n\n\n\n\n","category":"function"},{"location":"transform/#Equinoctial-Orbital-Elements-(eoe)","page":"Trasformations","title":"Equinoctial Orbital Elements (eoe)","text":"","category":"section"},{"location":"transform/","page":"Trasformations","title":"Trasformations","text":"convert6_cart_to_equi\nconvert6_equi_to_cart\nconvert6_coe_to_equi\nconvert6_equi_to_coe\n\n∂convert6_coe_to_equi","category":"page"},{"location":"transform/#AstroRepresentations.convert6_cart_to_equi","page":"Trasformations","title":"AstroRepresentations.convert6_cart_to_equi","text":"convert6_cart_to_equi(sv::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert cartesian state vector to Equinoctial keplerian elements.\n\nInputs\n\nsv – state vector – L, T\nμ – center gravitational parameter  – L³/T²\n\nOutput\n\nEquinoctial representation of the state as a SVector{6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications, 2013.\nCefola - DOI: 10.2514/6.1972-937\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.convert6_equi_to_cart","page":"Trasformations","title":"AstroRepresentations.convert6_equi_to_cart","text":"convert6_equi_to_cart(equi::AbstractVector{<:Number}, μ::Number, [args]...)\n\nConvert equinoctial state elements to cartesian state.\n\nInputs\n\nequi – Equinoctial elements – L, rad\nμ – center gravitational parameter  – L³/T²\n\nOutput\n\nCartesian representation of the state as a SVector{6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications, 2013.\nCefola - DOI: 10.2514/6.1972-937\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.convert6_coe_to_equi","page":"Trasformations","title":"AstroRepresentations.convert6_coe_to_equi","text":"convert6_coe_to_equi(coe::AbstractVector{<:Number}, [args]...)\n\nConvert classical orbital elements state vector to Equinoctial keplerian elements.\n\nInputs\n\ncoe – Keplerian elements – L, rad\n\nOutput\n\nEquinoctial representation of the state as a SVector{6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications, 2013.\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.convert6_equi_to_coe","page":"Trasformations","title":"AstroRepresentations.convert6_equi_to_coe","text":"convert6_equi_to_coe(coe::AbstractVector{<:Number}, [args]...)\n\nConvert Equinoctial keplerian elements to classical orbital elements state vector .\n\nInputs\n\nequi – Equinoctial elements – L, rad\n\nOutput\n\nKeplerian representation of the state as a SVector{6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications, 2013.\n\n\n\n\n\n","category":"function"},{"location":"transform/#AstroRepresentations.∂convert6_coe_to_equi","page":"Trasformations","title":"AstroRepresentations.∂convert6_coe_to_equi","text":"∂convert6_coe_to_equi(coe::AbstractVector{<:Number}, [args]...)\n\nConvert classical orbital elements state vector to Equinoctial keplerian elements.  Compute also the full jacobian of the equinoctial elemenents wrt the classical orbital elements.\n\nInputs\n\ncoe – Keplerian elements – L, rad\n\nOutput\n\nEquinoctial representation of the state as a SVector{6} and its jacobian as  a SMatrix{6, 6}.\n\nReferences\n\nVallado, David A. - Fundamentals of astrodynamics and applications, 2013.\n\n\n\n\n\n","category":"function"},{"location":"#AstroRepresentations.jl","page":"AstroRepresentations.jl","title":"AstroRepresentations.jl","text":"","category":"section"},{"location":"","page":"AstroRepresentations.jl","title":"AstroRepresentations.jl","text":"Astrodynamical states representations and transformations.","category":"page"}]
}
